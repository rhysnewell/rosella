use std::cell::RefCell;
use std::ffi::{c_int, c_void};
use std::ptr::NonNull;

use crate::learned::gemm::Layout;
use objc2::AnyThread;
use objc2::rc::{Retained, autoreleasepool};
use objc2::runtime::ProtocolObject;
use objc2_metal::{
    MTLCommandBuffer, MTLCommandBufferStatus, MTLCommandQueue, MTLCreateSystemDefaultDevice,
    MTLDevice, MTLResourceOptions,
};
use objc2_metal_performance_shaders::{
    MPSDataType, MPSMatrix, MPSMatrixDescriptor, MPSMatrixMultiplication,
};

// Each call maps its tables for the GPU, which costs about 1 ms. Only a product at least this
// long on every side ran faster there than on the matrix unit on an M1 Pro, waiting on it.
const SMALLEST: usize = 1024;

unsafe extern "C" {
    safe fn getpagesize() -> c_int;
}

struct Metal {
    device: Retained<ProtocolObject<dyn MTLDevice>>,
    queue: Retained<ProtocolObject<dyn MTLCommandQueue>>,
}

// A Metal buffer cannot cross threads, so the device lives on the thread that trains.
thread_local! {
    static METAL: RefCell<Option<Metal>> = const { RefCell::new(None) };
}

pub fn open() -> Option<String> {
    let device = MTLCreateSystemDefaultDevice()?;
    let queue = device.newCommandQueue()?;
    let name = device.name().to_string();
    METAL.set(Some(Metal { device, queue }));
    Some(name)
}

pub fn close() {
    METAL.set(None);
}

type End<P> = (P, Layout);

pub struct Running(Retained<ProtocolObject<dyn MTLCommandBuffer>>);

impl Drop for Running {
    fn drop(&mut self) {
        self.0.waitUntilCompleted();
        assert_eq!(
            self.0.status(),
            MTLCommandBufferStatus::Completed,
            "the GPU failed a product mid-run"
        );
    }
}

// A caller that goes on while the GPU works does not pay its fixed cost, so any size will do.
pub fn start(
    alpha: f32,
    ends: (End<*const f32>, End<*const f32>, End<*mut f32>),
    beta: f32,
    shape: (usize, usize, usize),
    any_size: bool,
) -> Option<Running> {
    let (m, n, k) = shape;
    if !any_size && m.min(n).min(k) < SMALLEST {
        return None;
    }
    METAL.with_borrow(|metal| metal.as_ref()?.start(alpha, ends, beta, shape))
}

impl Metal {
    // The buffer wraps the table where it lies. Only whole pages can be wrapped, and every page
    // from the one holding the first entry to the one holding the last is already mapped.
    fn matrix(&self, at: *const f32, layout: Layout) -> Option<Retained<MPSMatrix>> {
        let page = getpagesize() as usize;
        let start = at as usize;
        let base = start & !(page - 1);
        let end = (start + layout.elements() * size_of::<f32>()).next_multiple_of(page);
        let buffer = unsafe {
            self.device
                .newBufferWithBytesNoCopy_length_options_deallocator(
                    NonNull::new(base as *mut c_void)?,
                    end - base,
                    MTLResourceOptions::StorageModeShared,
                    None,
                )
        }?;
        let descriptor = unsafe {
            MPSMatrixDescriptor::matrixDescriptorWithRows_columns_rowBytes_dataType(
                layout.rows,
                layout.columns,
                layout.stride * size_of::<f32>(),
                MPSDataType::Float32,
            )
        };
        Some(unsafe {
            MPSMatrix::initWithBuffer_offset_descriptor(
                MPSMatrix::alloc(),
                &buffer,
                start - base,
                &descriptor,
            )
        })
    }

    fn start(
        &self,
        alpha: f32,
        ((a, left), (b, right), (c, out)): (End<*const f32>, End<*const f32>, End<*mut f32>),
        beta: f32,
        (m, n, k): (usize, usize, usize),
    ) -> Option<Running> {
        autoreleasepool(|_| {
            let (left_matrix, right_matrix, out_matrix) = (
                self.matrix(a, left)?,
                self.matrix(b, right)?,
                self.matrix(c, out)?,
            );
            let commands = self.queue.commandBuffer()?;
            unsafe {
                let kernel = MPSMatrixMultiplication::initWithDevice_transposeLeft_transposeRight_resultRows_resultColumns_interiorColumns_alpha_beta(
                    MPSMatrixMultiplication::alloc(),
                    &self.device,
                    left.turned,
                    right.turned,
                    m,
                    n,
                    k,
                    f64::from(alpha),
                    f64::from(beta),
                );
                kernel.encodeToCommandBuffer_leftMatrix_rightMatrix_resultMatrix(
                    &commands,
                    &left_matrix,
                    &right_matrix,
                    &out_matrix,
                );
            }
            commands.commit();
            Some(Running(commands))
        })
    }
}
