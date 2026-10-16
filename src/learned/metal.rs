use std::cell::RefCell;
use std::collections::HashMap;
use std::ffi::{c_int, c_void};
use std::ptr::NonNull;

use crate::learned::gemm::Layout;
use objc2::rc::{Retained, autoreleasepool};
use objc2::runtime::ProtocolObject;
use objc2::{AnyThread, msg_send};
use objc2_metal::{
    MTLCommandBuffer, MTLCommandBufferStatus, MTLCommandQueue, MTLCreateSystemDefaultDevice,
    MTLDevice, MTLResourceOptions,
};
use objc2_metal_performance_shaders::{
    MPSDataType, MPSMatrix, MPSMatrixDescriptor, MPSMatrixMultiplication,
};

// Below this on any side MPS gave a different sum each run, and it did not beat the matrix unit
// on an M1 Pro once each call's 1 ms of table mapping was paid.
const SMALLEST: usize = 1024;

#[repr(C)]
#[derive(Default)]
struct Timebase {
    numer: u32,
    denom: u32,
}

unsafe extern "C" {
    safe fn getpagesize() -> c_int;
    safe fn mach_absolute_time() -> u64;
    fn mach_timebase_info(info: *mut Timebase) -> c_int;
}

// The clock Metal stamps a command buffer's GPU start and end with.
fn host_seconds() -> f64 {
    let mut base = Timebase::default();
    unsafe { mach_timebase_info(&mut base) };
    mach_absolute_time() as f64 * f64::from(base.numer) / f64::from(base.denom) * 1e-9
}

type Key = (usize, usize, bool, bool);

struct Metal {
    device: Retained<ProtocolObject<dyn MTLDevice>>,
    queue: Retained<ProtocolObject<dyn MTLCommandQueue>>,
    shares: RefCell<HashMap<Key, f64>>,
}

// A Metal buffer cannot cross threads, so the device lives on the thread that trains.
thread_local! {
    static METAL: RefCell<Option<Metal>> = const { RefCell::new(None) };
}

pub fn open() -> Option<String> {
    let device = MTLCreateSystemDefaultDevice()?;
    let queue = device.newCommandQueue()?;
    let name = device.name().to_string();
    METAL.set(Some(Metal {
        device,
        queue,
        shares: RefCell::new(HashMap::new()),
    }));
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
) -> Option<Running> {
    METAL.with_borrow(|metal| metal.as_ref()?.start(alpha, ends, beta, shape))
}

pub struct Split {
    running: Running,
    key: Key,
    rows: usize,
    total: usize,
    started: f64,
}

// The GPU takes the first rows while the caller runs the rest on the matrix unit. Both sum each
// cell in the same order, so where the cut falls never changes a bit of the product.
pub fn split(
    alpha: f32,
    ((a, left), right, (c, out)): (End<*const f32>, End<*const f32>, End<*mut f32>),
    beta: f32,
    (m, n, k): (usize, usize, usize),
) -> Option<Split> {
    if m.min(n).min(k) < SMALLEST {
        return None;
    }
    let key = (n, k, left.turned, right.1.turned);
    METAL.with_borrow(|metal| {
        let metal = metal.as_ref()?;
        let share = metal.shares.borrow().get(&key).copied().unwrap_or(0.5);
        let rows = ((share * m as f64).round() as usize).clamp(SMALLEST, m);
        let ends = ((a, left.first(rows)), right, (c, out.first(rows)));
        let running = metal.start(alpha, ends, beta, (rows, n, k))?;
        Some(Split {
            running,
            key,
            rows,
            total: m,
            started: host_seconds(),
        })
    })
}

impl Split {
    pub fn rows(&self) -> usize {
        self.rows
    }

    // Moves the cut by half the rows that would have levelled this call, since the forward and
    // backward calls under one key find the GPU queue in different states.
    pub fn finish(self) {
        let done = host_seconds();
        self.running.0.waitUntilCompleted();
        let commands = &*self.running.0;
        let (start, end): (f64, f64) = unsafe {
            (
                msg_send![commands, GPUStartTime],
                msg_send![commands, GPUEndTime],
            )
        };
        let gpu_row = (end - start).max(0.0) / self.rows as f64;
        let cpu_row = match self.total - self.rows {
            0 => gpu_row,
            rows => (done - self.started) / rows as f64,
        };
        let moved = (end - done) / (gpu_row + cpu_row).max(f64::MIN_POSITIVE);
        let share = ((self.rows as f64 - moved / 2.0) / self.total as f64).clamp(0.0, 1.0);
        METAL.with_borrow(|metal| {
            if let Some(metal) = metal {
                metal.shares.borrow_mut().insert(self.key, share);
            }
        });
    }
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
