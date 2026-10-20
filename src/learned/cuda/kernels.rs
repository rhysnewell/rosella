use std::sync::Arc;

use anyhow::Result;
use cudarc::cublas::{result::sgemm, sys::cublasOperation_t};
use cudarc::driver::{CudaFunction, CudaModule, LaunchConfig, PushKernelArg};

use crate::learned::cuda::Device;

pub(crate) const SOURCE: &str = include_str!("kernels.cu");
const THREADS: u32 = 256;

macro_rules! kernels {
    ($($name:ident($($arg:ident: $kind:ty),*);)*) => {
        pub(crate) struct Kernels {
            $($name: CudaFunction,)*
        }

        impl Kernels {
            pub(crate) fn load(module: &Arc<CudaModule>) -> Result<Self> {
                Ok(Self {
                    $($name: module.load_function(stringify!($name))?,)*
                })
            }
        }

        impl Device {
            $(
                #[allow(clippy::too_many_arguments)] // one argument per kernel parameter
                pub(crate) fn $name(&self, config: LaunchConfig, $($arg: $kind),*) -> Result<()> {
                    let mut launch = self.stream.launch_builder(&self.kernels.$name);
                    $(launch.arg(&$arg);)*
                    unsafe { launch.launch(config) }?;
                    Ok(())
                }
            )*
        }
    };
}

kernels! {
    block_sums(x: u64, y: u64, mean: u64, rows: i32, cols: i32, block: i32, mode: i32, partial: u64);
    fold(partial: u64, blocks: i32, cols: i32, divisor: f64, out64: u64, out32: u64);
    finish(partial: u64, blocks: i32, cols: i32, rows: i32, mean: u64, momentum: f32, epsilon: f64,
        running_mean: u64, running_variance: u64, inverse: u64, mean32: u64);
    normalise(table: u64, active: u64, mean: u64, inverse: u64, scale: u64, shift: u64, rows: i32,
        cols: i32, dropping: i32, bar: u64, key: u64, kept: f32);
    mask(gradient: u64, active: u64, count: i64, kept: f32);
    rewrite(gradient: u64, normed: u64, scale: u64, inverse: u64, shift: u64, scaled: u64,
        rows: i32, cols: i32);
    broadcast(table: u64, row: u64, rows: i32, cols: i32);
    unit(table: u64, lengths: u64, rows: i32, cols: i32);
    unit_backward(gradient: u64, output: u64, lengths: u64, rows: i32, cols: i32);
    join(features: u64, width: i32, split: i32, depth: u64, depth_cols: i32, lift: f32,
        joined: u64, rows: i32);
    scaled_copy(from: u64, to: u64, count: i64, by: f32);
    adamw(params: u64, grads: u64, first: u64, second: u64, count: i64, shrink: f32, step: f32,
        momentum: f32, keep: f32, correction: f32, epsilon: f32);
    per_positive(weights: u64, n: i32, group_start: u64, group_size: u64, members: u64,
        lanes: i32, row_loss: u64);
}

pub(crate) fn spread(count: usize) -> LaunchConfig {
    LaunchConfig {
        grid_dim: ((count as u32).div_ceil(THREADS).max(1), 1, 1),
        block_dim: (THREADS, 1, 1),
        shared_mem_bytes: 0,
    }
}

pub(crate) fn columns(cols: usize, blocks: usize) -> LaunchConfig {
    LaunchConfig {
        grid_dim: ((cols as u32).div_ceil(THREADS), blocks as u32, 1),
        ..spread(0)
    }
}

pub(crate) fn per_row(rows: usize) -> LaunchConfig {
    LaunchConfig {
        grid_dim: (rows as u32, 1, 1),
        ..spread(0)
    }
}

#[derive(Clone, Copy)]
pub(crate) struct Table {
    pub(crate) at: u64,
    pub(crate) rows: usize,
    pub(crate) cols: usize,
    pub(crate) stride: usize,
}

impl Table {
    pub(crate) fn new(at: u64, rows: usize, cols: usize) -> Self {
        Self {
            at,
            rows,
            cols,
            stride: cols,
        }
    }

    pub(crate) fn from_column(self, column: usize) -> Self {
        Self {
            at: self.at + (column * size_of::<f32>()) as u64,
            cols: self.cols - column,
            ..self
        }
    }
}

impl Device {
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
                b.at as *const f32,
                b.stride as i32,
                a.at as *const f32,
                a.stride as i32,
                &beta,
                c.at as *mut f32,
                c.stride as i32,
            )
        }?;
        Ok(())
    }
}
