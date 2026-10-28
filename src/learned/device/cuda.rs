use std::cell::RefCell;
use std::sync::Arc;

use anyhow::Result;
use cudarc::cublas::CudaBlas;
use cudarc::cublas::result::sgemm;
use cudarc::cublas::sys::cublasOperation_t;
use cudarc::driver::sys::CUdevice_attribute::CU_DEVICE_ATTRIBUTE_SINGLE_TO_DOUBLE_PRECISION_PERF_RATIO;
use cudarc::driver::{
    CudaContext, CudaFunction, CudaModule, CudaSlice, CudaStream, DevicePtr, DeviceRepr,
    LaunchConfig, PushKernelArg, result,
};
use cudarc::nvrtc::{CompileOptions, compile_ptx_with_opts};

use crate::learned::device::{Launch, Ptr, Table};

const SOURCE: &str = concat!(include_str!("kernels.cu"), "\n", include_str!("nearest.cu"));
// Data centre cards run doubles at a half or a third of the float rate. Gaming cards run them at a
// 32nd or a 64th, where the coverage distance in f64 would be slower than on the cores.
const SLOWEST_DOUBLES: i32 = 4;

unsafe impl DeviceRepr for Ptr {}

macro_rules! cuda_kernels {
    ($($name:ident($($arg:ident: $kind:ty),*);)*) => {
        struct Kernels {
            $($name: CudaFunction,)*
        }

        impl Kernels {
            fn load(module: &Arc<CudaModule>) -> Result<Self> {
                Ok(Self {
                    $($name: module.load_function(stringify!($name))?,)*
                })
            }
        }

        impl Device {
            $(
                #[allow(clippy::too_many_arguments)] // one argument per kernel parameter
                pub(crate) fn $name(&self, launch: Launch, $($arg: $kind),*) -> Result<()> {
                    let mut call = self.stream.launch_builder(&self.kernels.$name);
                    $(call.arg(&$arg);)*
                    unsafe { call.launch(config(launch)) }?;
                    Ok(())
                }
            )*
        }
    };
}

kernels!(cuda_kernels);

fn config(launch: Launch) -> LaunchConfig {
    LaunchConfig {
        grid_dim: (launch.groups.0, launch.groups.1, 1),
        block_dim: (launch.threads, 1, 1),
        shared_mem_bytes: 0,
    }
}

pub(crate) struct Device {
    stream: Arc<CudaStream>,
    blas: CudaBlas,
    kernels: Kernels,
    held: RefCell<Vec<CudaSlice<u8>>>,
    wide: bool,
}

impl Device {
    // Every library loads at run time, so a machine without CUDA runs the same binary on its cores.
    pub(crate) fn missing() -> Result<Option<String>> {
        let found = unsafe {
            [
                ("libcuda", cudarc::driver::sys::is_culib_present()),
                ("libcublas", cudarc::cublas::sys::is_culib_present()),
                ("libnvrtc", cudarc::nvrtc::sys::is_culib_present()),
            ]
        };
        if let Some((library, _)) = found.iter().find(|(_, here)| !here) {
            return Ok(Some(format!("{library} not found")));
        }
        Ok((CudaContext::device_count()? == 0).then(|| "no CUDA device".to_string()))
    }

    pub(crate) fn new() -> Result<(Self, String)> {
        let context = CudaContext::new(0)?;
        let (major, minor) = context.compute_capability()?;
        let wide = context.attribute(CU_DEVICE_ATTRIBUTE_SINGLE_TO_DOUBLE_PRECISION_PERF_RATIO)?
            <= SLOWEST_DOUBLES;
        let ptx = compile_ptx_with_opts(
            SOURCE,
            CompileOptions {
                options: vec![
                    format!("--gpu-architecture=compute_{major}{minor}"),
                    format!("-DNORM_EPSILON={:e}", crate::learned::norm::EPSILON),
                    format!("-DCOVERAGE={}", if wide { "double" } else { "float" }),
                ],
                ..Default::default()
            },
        )?;
        let module = context.load_module(ptx)?;
        let stream = context.default_stream();
        let blas = CudaBlas::new(stream.clone())?;
        let device = Self {
            kernels: Kernels::load(&module)?,
            stream,
            blas,
            held: RefCell::new(Vec::new()),
            wide,
        };
        Ok((device, context.name()?))
    }

    pub(crate) fn wide(&self) -> bool {
        self.wide
    }

    pub(crate) fn alloc(&self, bytes: usize) -> Result<Ptr> {
        let slice = self.stream.alloc_zeros::<u8>(bytes.max(1))?;
        let at = Ptr(slice.device_ptr(&self.stream).0);
        self.held.borrow_mut().push(slice);
        Ok(at)
    }

    pub(crate) fn mark(&self) -> usize {
        self.held.borrow().len()
    }

    pub(crate) fn release(&self, mark: usize) -> Result<()> {
        self.stream.synchronize()?;
        self.held.borrow_mut().truncate(mark);
        Ok(())
    }

    pub(crate) fn upload<T>(&self, at: Ptr, values: &[T]) -> Result<()> {
        unsafe { result::memcpy_htod_async(at.0, values, self.stream.cu_stream()) }?;
        Ok(())
    }

    pub(crate) fn download<T: Copy + Default>(&self, at: Ptr, len: usize) -> Result<Vec<T>> {
        let mut values = vec![T::default(); len];
        unsafe { result::memcpy_dtoh_async(&mut values, at.0, self.stream.cu_stream()) }?;
        self.stream.synchronize()?;
        Ok(values)
    }

    pub(crate) fn sums(&self, at: Ptr, len: usize) -> Result<Vec<f64>> {
        self.download(at, len)
    }

    // cuBLAS reads a row-major table as its turn. So the product is taken turned, operands swapped.
    pub(crate) fn product(
        &self,
        alpha: f32,
        (a, a_turned): (Table, bool),
        (b, b_turned): (Table, bool),
        beta: f32,
        c: Table,
    ) -> Result<()> {
        let operation = |turned: bool| match turned {
            true => cublasOperation_t::CUBLAS_OP_T,
            false => cublasOperation_t::CUBLAS_OP_N,
        };
        let inner = if a_turned { a.rows } else { a.cols };
        unsafe {
            sgemm(
                *self.blas.handle(),
                operation(b_turned),
                operation(a_turned),
                c.cols as i32,
                c.rows as i32,
                inner as i32,
                &alpha,
                b.at.0 as *const f32,
                b.stride as i32,
                a.at.0 as *const f32,
                a.stride as i32,
                &beta,
                c.at.0 as *mut f32,
                c.stride as i32,
            )
        }?;
        Ok(())
    }
}
