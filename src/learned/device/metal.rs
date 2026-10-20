use std::cell::RefCell;
use std::collections::HashMap;
use std::ffi::c_void;
use std::ptr::NonNull;

use anyhow::{Result, anyhow, bail};
use objc2::AnyThread;
use objc2::rc::{Retained, autoreleasepool};
use objc2::runtime::ProtocolObject;
use objc2_foundation::NSString;
use objc2_metal::{
    MTLBuffer, MTLCommandBuffer, MTLCommandBufferStatus, MTLCommandEncoder, MTLCommandQueue,
    MTLCompileOptions, MTLComputeCommandEncoder, MTLComputePipelineState,
    MTLCreateSystemDefaultDevice, MTLDevice, MTLLibrary, MTLResourceOptions, MTLSize,
};
use objc2_metal_performance_shaders::{
    MPSDataType, MPSMatrix, MPSMatrixDescriptor, MPSMatrixMultiplication,
};

use crate::learned::device::{Launch, Ptr, Table, spread};

const SOURCE: &str = include_str!("kernels.metal");
// A handle is the buffer's index above this bit and a byte offset below, so the shared code adds
// offsets to it as it does to a CUDA pointer.
const SHIFT: u32 = 40;
const SKINNY: usize = 16;

type Buffer = Retained<ProtocolObject<dyn MTLBuffer>>;
type Pipeline = Retained<ProtocolObject<dyn MTLComputePipelineState>>;
type Encoder = Retained<ProtocolObject<dyn MTLComputeCommandEncoder>>;
type Commands = Retained<ProtocolObject<dyn MTLCommandBuffer>>;
type Key = (bool, bool, usize, usize, usize, u32, u32);

trait Bind {
    fn bind(&self, device: &Device, encoder: &Encoder, index: usize);
}

impl Bind for Ptr {
    fn bind(&self, device: &Device, encoder: &Encoder, index: usize) {
        let (buffer, offset) = device.resolve(*self);
        unsafe { encoder.setBuffer_offset_atIndex(Some(&buffer), offset, index) };
    }
}

macro_rules! scalar {
    ($($kind:ty),*) => {$(
        impl Bind for $kind {
            fn bind(&self, _: &Device, encoder: &Encoder, index: usize) {
                let bytes = NonNull::from(self).cast::<c_void>();
                unsafe { encoder.setBytes_length_atIndex(bytes, size_of::<$kind>(), index) };
            }
        }
    )*};
}

scalar!(i32, u64, f32);

macro_rules! metal_kernels {
    ($($name:ident($($arg:ident: $kind:ty),*);)*) => {
        struct Kernels {
            $($name: Pipeline,)*
            dot: Pipeline,
        }

        impl Kernels {
            fn load(
                device: &ProtocolObject<dyn MTLDevice>,
                library: &ProtocolObject<dyn MTLLibrary>,
            ) -> Result<Self> {
                Ok(Self {
                    $($name: pipeline(device, library, stringify!($name))?,)*
                    dot: pipeline(device, library, "dot")?,
                })
            }
        }

        impl Device {
            $(
                #[allow(clippy::too_many_arguments)] // one argument per kernel parameter
                pub(crate) fn $name(&self, launch: Launch, $($arg: $kind),*) -> Result<()> {
                    self.dispatch(&self.kernels.$name, launch, &[$(&$arg as &dyn Bind),*])
                }
            )*
        }
    };
}

kernels!(metal_kernels);

fn pipeline(
    device: &ProtocolObject<dyn MTLDevice>,
    library: &ProtocolObject<dyn MTLLibrary>,
    name: &str,
) -> Result<Pipeline> {
    let function = library
        .newFunctionWithName(&NSString::from_str(name))
        .ok_or_else(|| anyhow!("the Metal library has no {name}"))?;
    device
        .newComputePipelineStateWithFunction_error(&function)
        .map_err(|error| anyhow!("{name}: {error}"))
}

#[derive(Default)]
struct Pending {
    commands: Option<Commands>,
    encoder: Option<Encoder>,
}

pub(crate) struct Device {
    device: Retained<ProtocolObject<dyn MTLDevice>>,
    queue: Retained<ProtocolObject<dyn MTLCommandQueue>>,
    kernels: Kernels,
    buffers: RefCell<Vec<Buffer>>,
    pending: RefCell<Pending>,
    products: RefCell<HashMap<Key, Retained<MPSMatrixMultiplication>>>,
}

impl Device {
    pub(crate) fn missing() -> Result<Option<String>> {
        Ok(MTLCreateSystemDefaultDevice()
            .is_none()
            .then(|| "no Metal device".to_string()))
    }

    // Fast math assumes no infinities, and the loss masks a row's own entry with one.
    pub(crate) fn new() -> Result<(Self, String)> {
        let device = MTLCreateSystemDefaultDevice().ok_or_else(|| anyhow!("no Metal device"))?;
        let queue = device
            .newCommandQueue()
            .ok_or_else(|| anyhow!("no Metal command queue"))?;
        let options = MTLCompileOptions::new();
        #[allow(deprecated)] // mathMode needs macOS 15, and this switch works on every release
        options.setFastMathEnabled(false);
        let library = device
            .newLibraryWithSource_options_error(&NSString::from_str(SOURCE), Some(&options))
            .map_err(|error| anyhow!("the Metal kernels did not compile: {error}"))?;
        let kernels = Kernels::load(&device, &library)?;
        let name = device.name().to_string();
        let null = device
            .newBufferWithLength_options(16, MTLResourceOptions::StorageModeShared)
            .ok_or_else(|| anyhow!("no Metal buffer"))?;
        Ok((
            Self {
                device,
                queue,
                kernels,
                buffers: RefCell::new(vec![null]),
                pending: RefCell::new(Pending::default()),
                products: RefCell::new(HashMap::new()),
            },
            name,
        ))
    }

    fn resolve(&self, at: Ptr) -> (Buffer, usize) {
        let index = (at.0 >> SHIFT) as usize;
        let offset = (at.0 & ((1 << SHIFT) - 1)) as usize;
        (self.buffers.borrow()[index].clone(), offset)
    }

    fn host(&self, at: Ptr) -> *mut u8 {
        let (buffer, offset) = self.resolve(at);
        unsafe { buffer.contents().cast::<u8>().as_ptr().add(offset) }
    }

    pub(crate) fn alloc(&self, bytes: usize) -> Result<Ptr> {
        let buffer = self
            .device
            .newBufferWithLength_options(bytes.max(16), MTLResourceOptions::StorageModeShared)
            .ok_or_else(|| anyhow!("Metal could not hold {bytes} bytes"))?;
        unsafe { std::ptr::write_bytes(buffer.contents().cast::<u8>().as_ptr(), 0, bytes) };
        let mut buffers = self.buffers.borrow_mut();
        buffers.push(buffer);
        Ok(Ptr(((buffers.len() - 1) as u64) << SHIFT))
    }

    // Memory is shared, so a copy is a plain one once the device has finished what it was given.
    pub(crate) fn upload<T: Copy>(&self, at: Ptr, values: &[T]) -> Result<()> {
        self.flush()?;
        unsafe {
            std::ptr::copy_nonoverlapping(values.as_ptr(), self.host(at).cast::<T>(), values.len())
        };
        Ok(())
    }

    pub(crate) fn download<T: Copy + Default>(&self, at: Ptr, len: usize) -> Result<Vec<T>> {
        self.flush()?;
        let mut values = vec![T::default(); len];
        unsafe {
            std::ptr::copy_nonoverlapping(self.host(at).cast::<T>(), values.as_mut_ptr(), len)
        };
        Ok(values)
    }

    fn commands(&self) -> Result<Commands> {
        let mut pending = self.pending.borrow_mut();
        if pending.commands.is_none() {
            let commands = self
                .queue
                .commandBuffer()
                .ok_or_else(|| anyhow!("no Metal command buffer"))?;
            pending.commands = Some(commands);
        }
        Ok(pending.commands.clone().expect("opened above"))
    }

    fn encoder(&self) -> Result<Encoder> {
        let commands = self.commands()?;
        let mut pending = self.pending.borrow_mut();
        if pending.encoder.is_none() {
            let encoder = commands
                .computeCommandEncoder()
                .ok_or_else(|| anyhow!("no Metal compute encoder"))?;
            pending.encoder = Some(encoder);
        }
        Ok(pending.encoder.clone().expect("opened above"))
    }

    fn close_encoder(&self) {
        if let Some(encoder) = self.pending.borrow_mut().encoder.take() {
            encoder.endEncoding();
        }
    }

    fn flush(&self) -> Result<()> {
        self.close_encoder();
        let Some(commands) = self.pending.borrow_mut().commands.take() else {
            return Ok(());
        };
        commands.commit();
        commands.waitUntilCompleted();
        if commands.status() != MTLCommandBufferStatus::Completed {
            bail!("the GPU failed a training step");
        }
        Ok(())
    }

    fn dispatch(&self, pipeline: &Pipeline, launch: Launch, args: &[&dyn Bind]) -> Result<()> {
        let encoder = self.encoder()?;
        encoder.setComputePipelineState(pipeline);
        for (index, arg) in args.iter().enumerate() {
            arg.bind(self, &encoder, index);
        }
        let size = |width, height| MTLSize {
            width,
            height,
            depth: 1,
        };
        encoder.dispatchThreadgroups_threadsPerThreadgroup(
            size(launch.groups.0 as usize, launch.groups.1 as usize),
            size(launch.threads as usize, 1),
        );
        Ok(())
    }

    fn matrix(&self, table: Table) -> Retained<MPSMatrix> {
        let (buffer, offset) = self.resolve(table.at);
        let descriptor = unsafe {
            MPSMatrixDescriptor::matrixDescriptorWithRows_columns_rowBytes_dataType(
                table.rows,
                table.cols,
                table.stride * size_of::<f32>(),
                MPSDataType::Float32,
            )
        };
        unsafe {
            MPSMatrix::initWithBuffer_offset_descriptor(
                MPSMatrix::alloc(),
                &buffer,
                offset,
                &descriptor,
            )
        }
    }

    pub(crate) fn product(
        &self,
        alpha: f32,
        (a, a_turned): (Table, bool),
        (b, b_turned): (Table, bool),
        beta: f32,
        c: Table,
    ) -> Result<()> {
        let inner = if a_turned { a.rows } else { a.cols };
        if c.rows.min(c.cols).min(inner) < SKINNY {
            let strides = |table: Table, turned: bool| match turned {
                true => (1, table.stride as i32),
                false => (table.stride as i32, 1),
            };
            let ((a_rows, a_cols), (b_rows, b_cols)) = (strides(a, a_turned), strides(b, b_turned));
            let (m, n, k) = (c.rows as i32, c.cols as i32, inner as i32);
            return self.dispatch(
                &self.kernels.dot,
                spread(c.rows * c.cols),
                &[
                    &a.at,
                    &a_rows,
                    &a_cols,
                    &b.at,
                    &b_rows,
                    &b_cols,
                    &c.at,
                    &(c.stride as i32),
                    &m,
                    &n,
                    &k,
                    &alpha,
                    &beta,
                ],
            );
        }
        self.close_encoder();
        let commands = self.commands()?;
        let key = (
            a_turned,
            b_turned,
            c.rows,
            c.cols,
            inner,
            alpha.to_bits(),
            beta.to_bits(),
        );
        let kernel = self
            .products
            .borrow_mut()
            .entry(key)
            .or_insert_with(|| unsafe {
                MPSMatrixMultiplication::initWithDevice_transposeLeft_transposeRight_resultRows_resultColumns_interiorColumns_alpha_beta(
                    MPSMatrixMultiplication::alloc(),
                    &self.device,
                    a_turned,
                    b_turned,
                    c.rows,
                    c.cols,
                    inner,
                    f64::from(alpha),
                    f64::from(beta),
                )
            })
            .clone();
        autoreleasepool(|_| unsafe {
            kernel.encodeToCommandBuffer_leftMatrix_rightMatrix_resultMatrix(
                &commands,
                &self.matrix(a),
                &self.matrix(b),
                &self.matrix(c),
            );
        });
        Ok(())
    }
}

impl Drop for Device {
    fn drop(&mut self) {
        let _ = self.flush();
    }
}
