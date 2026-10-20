mod kernels;
mod loss;
mod net;
mod trainer;

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use anyhow::Result;
use cudarc::cublas::CudaBlas;
use cudarc::driver::{CudaContext, CudaSlice, CudaStream, DevicePtr, DeviceRepr, ValidAsZeroBits};
use cudarc::nvrtc::{CompileOptions, compile_ptx_with_opts};

use kernels::{Kernels, SOURCE};
pub use trainer::Trainer;

pub(crate) struct Device {
    stream: Arc<CudaStream>,
    blas: CudaBlas,
    kernels: Kernels,
}

// The trainer reads the device off the thread that opened it, as the Metal path does.
thread_local! {
    static DEVICE: RefCell<Option<Rc<Device>>> = const { RefCell::new(None) };
}

pub enum Found {
    Device(String),
    Missing(String),
}

// Every library loads at run time, so a machine without CUDA runs the same binary on its cores.
fn missing() -> Result<Option<String>> {
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

// No fused multiply-add, so each kernel rounds as the host passes do.
fn device() -> Result<(Device, String)> {
    let context = CudaContext::new(0)?;
    let (major, minor) = context.compute_capability()?;
    let ptx = compile_ptx_with_opts(
        SOURCE,
        CompileOptions {
            fmad: Some(false),
            options: vec![format!("--gpu-architecture=compute_{major}{minor}")],
            ..Default::default()
        },
    )?;
    let module = context.load_module(ptx)?;
    let stream = context.default_stream();
    let blas = CudaBlas::new(stream.clone())?;
    let name = context.name()?;
    Ok((
        Device {
            kernels: Kernels::load(&module)?,
            stream,
            blas,
        },
        name,
    ))
}

pub fn open() -> Result<Found> {
    if let Some(why) = missing()? {
        return Ok(Found::Missing(why));
    }
    let (device, name) = device()?;
    DEVICE.set(Some(Rc::new(device)));
    Ok(Found::Device(name))
}

pub fn close() {
    DEVICE.set(None);
}

pub(crate) fn current() -> Option<Rc<Device>> {
    DEVICE.with_borrow(Clone::clone)
}

impl Device {
    pub(crate) fn zeros<T: DeviceRepr + ValidAsZeroBits>(
        &self,
        len: usize,
    ) -> Result<CudaSlice<T>> {
        Ok(self.stream.alloc_zeros(len.max(1))?)
    }

    pub(crate) fn upload<T: DeviceRepr>(&self, from: &[T], to: &mut CudaSlice<T>) -> Result<()> {
        let mut view = to.slice_mut(..from.len());
        Ok(self.stream.memcpy_htod(from, &mut view)?)
    }

    pub(crate) fn download<T: DeviceRepr + Default + Clone>(
        &self,
        from: &CudaSlice<T>,
        len: usize,
    ) -> Result<Vec<T>> {
        let mut to = vec![T::default(); len];
        self.stream.memcpy_dtoh(&from.slice(..len), &mut to)?;
        Ok(to)
    }

    pub(crate) fn at<T>(&self, slice: &CudaSlice<T>) -> u64 {
        slice.device_ptr(&self.stream).0
    }
}
