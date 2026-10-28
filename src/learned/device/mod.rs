macro_rules! kernels {
    ($apply:ident) => {
        $apply! {
            block_sums(x: Ptr, y: Ptr, mean: Ptr, rows: i32, cols: i32, block: i32, mode: i32,
                partial: Ptr);
            fold(partial: Ptr, blocks: i32, cols: i32, divisor: f32, out: Ptr);
            finish(partial: Ptr, blocks: i32, cols: i32, rows: i32, mean: Ptr, momentum: f32,
                epsilon: f32, running_mean: Ptr, running_variance: Ptr, inverse: Ptr);
            normalise(table: Ptr, active: Ptr, mean: Ptr, inverse: Ptr, scale: Ptr, shift: Ptr,
                rows: i32, cols: i32, dropping: i32, bar: u64, key: u64, kept: f32);
            mask(gradient: Ptr, active: Ptr, count: i32, kept: f32);
            rewrite(gradient: Ptr, normed: Ptr, scale: Ptr, inverse: Ptr, shift: Ptr, scaled: Ptr,
                rows: i32, cols: i32);
            broadcast(table: Ptr, row: Ptr, rows: i32, cols: i32);
            unit(table: Ptr, lengths: Ptr, rows: i32, cols: i32);
            unit_backward(gradient: Ptr, output: Ptr, lengths: Ptr, rows: i32, cols: i32);
            join(features: Ptr, width: i32, split: i32, depth: Ptr, depth_cols: i32, lift: f32,
                joined: Ptr, rows: i32);
            scaled_copy(from: Ptr, to: Ptr, count: i32, by: f32);
            adamw(params: Ptr, grads: Ptr, first: Ptr, second: Ptr, count: i32, shrink: f32,
                step: f32, momentum: f32, keep: f32, correction: f32, epsilon: f32);
            per_positive(weights: Ptr, n: i32, group_start: Ptr, group_size: Ptr, members: Ptr,
                peak: f32, scale: f32, row_loss: Ptr);
            contrastive(weights: Ptr, n: i32, group_start: Ptr, group_size: Ptr, members: Ptr,
                peak: f32, scale: f32, row_loss: Ptr);
            nearest(dots: Ptr, rows: i32, cols: i32, first: i32, selves: Ptr, query_variance: Ptr,
                base_variance: Ptr, k: i32, distances: Ptr, indices: Ptr, query_rows: Ptr,
                moments: Ptr, presence: Ptr, stride: i32, samples: i32, terms: Ptr);
        }
    };
}

#[cfg(not(target_os = "macos"))]
mod cuda;
mod knn;
mod loss;
#[cfg(target_os = "macos")]
mod metal;
mod net;
mod trainer;

use std::cell::RefCell;
use std::rc::Rc;

use anyhow::Result;

#[cfg(not(target_os = "macos"))]
use cuda::Device;
pub use knn::{Coverage, NOBODY, Rows, nearest};
#[cfg(target_os = "macos")]
use metal::Device;
pub use trainer::{Single, Trainer, embed};

const THREADS: u32 = 256;

#[derive(Clone, Copy, Default)]
#[repr(transparent)]
pub(crate) struct Ptr(u64);

impl Ptr {
    fn after(self, entries: usize) -> Self {
        Self(self.0 + (entries * size_of::<f32>()) as u64)
    }
}

#[derive(Clone, Copy)]
pub(crate) struct Launch {
    groups: (u32, u32),
    threads: u32,
}

fn spread(count: usize) -> Launch {
    Launch {
        groups: ((count as u32).div_ceil(THREADS).max(1), 1),
        threads: THREADS,
    }
}

fn columns(cols: usize, blocks: usize) -> Launch {
    Launch {
        groups: ((cols as u32).div_ceil(THREADS), blocks as u32),
        threads: THREADS,
    }
}

fn per_row(rows: usize) -> Launch {
    Launch {
        groups: (rows as u32, 1),
        threads: THREADS,
    }
}

#[derive(Clone, Copy)]
pub(crate) struct Table {
    at: Ptr,
    rows: usize,
    cols: usize,
    stride: usize,
}

impl Table {
    fn new(at: Ptr, rows: usize, cols: usize) -> Self {
        Self {
            at,
            rows,
            cols,
            stride: cols,
        }
    }

    fn cols_from(self, column: usize) -> Self {
        Self {
            at: self.at.after(column),
            cols: self.cols - column,
            ..self
        }
    }
}

pub enum Found {
    Device(String),
    Missing(String),
}

// The trainer reads the device off the thread that opened it.
thread_local! {
    static DEVICE: RefCell<Option<Rc<Device>>> = const { RefCell::new(None) };
}

pub fn open() -> Result<Found> {
    if let Some(why) = Device::missing()? {
        return Ok(Found::Missing(why));
    }
    let (device, name) = Device::new()?;
    DEVICE.set(Some(Rc::new(device)));
    Ok(Found::Device(name))
}

pub fn close() {
    DEVICE.set(None);
}

pub fn is_open() -> bool {
    DEVICE.with_borrow(Option::is_some)
}

fn current() -> Option<Rc<Device>> {
    DEVICE.with_borrow(Clone::clone)
}
