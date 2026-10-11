use ndarray::{Array1, Array2, ArrayView2, s};
use rand::{Rng, rngs::StdRng};
use rayon::prelude::*;

use crate::learned::gemm::product;

pub const HIDDEN: usize = 512;
pub const WIDTH: usize = 64;
const EPSILON: f64 = 1e-5;
const MOMENTUM: f32 = 0.1;
// Fixed blocks keep every column sum in one order whatever the thread count.
const BLOCK: usize = 256;

#[derive(Clone)]
pub struct Params {
    first: Array2<f32>,
    first_scale: Array1<f32>,
    first_shift: Array1<f32>,
    second: Array2<f32>,
    second_scale: Array1<f32>,
    second_shift: Array1<f32>,
    last: Array2<f32>,
    last_bias: Array1<f32>,
}

impl Params {
    // The uniform bound PyTorch gives a linear layer by default.
    fn new(input: usize, rng: &mut StdRng) -> Self {
        let mut uniform = |rows: usize, columns: usize| {
            let bound = 1.0 / (columns as f32).sqrt();
            Array2::from_shape_simple_fn((rows, columns), || rng.random_range(-bound..bound))
        };
        let first = uniform(HIDDEN, input);
        let second = uniform(HIDDEN, HIDDEN);
        let last = uniform(WIDTH, HIDDEN);
        let last_bias = uniform(1, HIDDEN).slice(s![0, ..WIDTH]).to_owned();
        Self {
            first,
            first_scale: Array1::ones(HIDDEN),
            first_shift: Array1::zeros(HIDDEN),
            second,
            second_scale: Array1::ones(HIDDEN),
            second_shift: Array1::zeros(HIDDEN),
            last,
            last_bias,
        }
    }

    pub fn zeros_like(&self) -> Self {
        Self {
            first: Array2::zeros(self.first.dim()),
            first_scale: Array1::zeros(HIDDEN),
            first_shift: Array1::zeros(HIDDEN),
            second: Array2::zeros(self.second.dim()),
            second_scale: Array1::zeros(HIDDEN),
            second_shift: Array1::zeros(HIDDEN),
            last: Array2::zeros(self.last.dim()),
            last_bias: Array1::zeros(WIDTH),
        }
    }

    pub fn slices(&self) -> [&[f32]; 8] {
        [
            flat(&self.first),
            flat(&self.first_scale),
            flat(&self.first_shift),
            flat(&self.second),
            flat(&self.second_scale),
            flat(&self.second_shift),
            flat(&self.last),
            flat(&self.last_bias),
        ]
    }

    pub fn slices_mut(&mut self) -> [&mut [f32]; 8] {
        [
            flat_mut(&mut self.first),
            flat_mut(&mut self.first_scale),
            flat_mut(&mut self.first_shift),
            flat_mut(&mut self.second),
            flat_mut(&mut self.second_scale),
            flat_mut(&mut self.second_shift),
            flat_mut(&mut self.last),
            flat_mut(&mut self.last_bias),
        ]
    }
}

fn flat<D: ndarray::Dimension>(values: &ndarray::Array<f32, D>) -> &[f32] {
    values.as_slice().expect("parameters are contiguous")
}

fn flat_mut<D: ndarray::Dimension>(values: &mut ndarray::Array<f32, D>) -> &mut [f32] {
    values.as_slice_mut().expect("parameters are contiguous")
}

struct Running {
    mean: Array1<f32>,
    variance: Array1<f32>,
}

pub struct Net {
    pub params: Params,
    running: [Running; 2],
}

pub struct Pass {
    input: Array2<f32>,
    normed: [Array2<f32>; 2],
    inverse: [Array1<f32>; 2],
    active: [Array2<f32>; 2],
    pub output: Array2<f32>,
    lengths: Vec<f32>,
}

fn column_sums(rows: usize, columns: usize, add: impl Fn(usize, &mut [f64]) + Sync) -> Vec<f64> {
    (0..rows.div_ceil(BLOCK))
        .into_par_iter()
        .map(|block| {
            let mut sums = vec![0.0; columns];
            for row in block * BLOCK..((block + 1) * BLOCK).min(rows) {
                add(row, &mut sums);
            }
            sums
        })
        .collect::<Vec<_>>()
        .into_iter()
        .fold(vec![0.0; columns], |mut total, sums| {
            total
                .iter_mut()
                .zip(sums)
                .for_each(|(total, sum)| *total += sum);
            total
        })
}

fn row(table: &Array2<f32>, at: usize) -> &[f32] {
    let columns = table.ncols();
    &table.as_slice().expect("tables are contiguous")[at * columns..(at + 1) * columns]
}

fn rows_mut(table: &mut Array2<f32>) -> rayon::slice::ChunksMut<'_, f32> {
    let columns = table.ncols();
    table
        .as_slice_mut()
        .expect("tables are contiguous")
        .par_chunks_mut(columns)
}

// Batch statistics normalise in place and leave the running pair for embedding (Ioffe & Szegedy
// 2015), the running variance unbiased as PyTorch keeps it.
fn normalise(table: &mut Array2<f32>, running: &mut Running) -> Array1<f32> {
    let (rows, columns) = table.dim();
    let mean = column_sums(rows, columns, |at, sums| {
        sums.iter_mut()
            .zip(row(table, at))
            .for_each(|(sum, value)| *sum += f64::from(*value));
    })
    .into_iter()
    .map(|sum| sum / rows as f64)
    .collect::<Vec<_>>();
    let variance = column_sums(rows, columns, |at, sums| {
        sums.iter_mut()
            .zip(row(table, at))
            .zip(&mean)
            .for_each(|((sum, value), mean)| *sum += (f64::from(*value) - mean).powi(2));
    })
    .into_iter()
    .map(|sum| sum / rows as f64)
    .collect::<Vec<_>>();
    let unbiased = rows as f64 / (rows.max(2) - 1) as f64;
    for at in 0..columns {
        running.mean[at] = (1.0 - MOMENTUM) * running.mean[at] + MOMENTUM * mean[at] as f32;
        running.variance[at] =
            (1.0 - MOMENTUM) * running.variance[at] + MOMENTUM * (variance[at] * unbiased) as f32;
    }
    let inverse = variance
        .iter()
        .map(|variance| (1.0 / (variance + EPSILON).sqrt()) as f32)
        .collect::<Array1<f32>>();
    let mean = mean.iter().map(|mean| *mean as f32).collect::<Vec<_>>();
    rows_mut(table).for_each(|values| {
        for ((value, mean), inverse) in values.iter_mut().zip(&mean).zip(&inverse) {
            *value = (*value - mean) * inverse;
        }
    });
    inverse
}

fn activate(normed: &Array2<f32>, scale: &Array1<f32>, shift: &Array1<f32>) -> Array2<f32> {
    let mut active = normed.clone();
    rows_mut(&mut active).for_each(|values| {
        for ((value, scale), shift) in values.iter_mut().zip(scale).zip(shift) {
            *value = (*value * scale + shift).max(0.0);
        }
    });
    active
}

fn unit(table: &mut Array2<f32>) -> Vec<f32> {
    let columns = table.ncols();
    table
        .as_slice_mut()
        .expect("tables are contiguous")
        .par_chunks_mut(columns)
        .map(|values| {
            let length = values.iter().map(|value| value * value).sum::<f32>().sqrt();
            let length = length.max(1e-12);
            values.iter_mut().for_each(|value| *value /= length);
            length
        })
        .collect()
}

impl Net {
    pub fn new(input: usize, rng: &mut StdRng) -> Self {
        let running = || Running {
            mean: Array1::zeros(HIDDEN),
            variance: Array1::ones(HIDDEN),
        };
        Self {
            params: Params::new(input, rng),
            running: [running(), running()],
        }
    }

    pub fn forward(&mut self, input: Array2<f32>) -> Pass {
        let rows = input.nrows();
        let mut first = Array2::zeros((rows, HIDDEN));
        product(
            1.0,
            input.view(),
            self.params.first.t(),
            0.0,
            first.view_mut(),
        );
        let first_inverse = normalise(&mut first, &mut self.running[0]);
        let first_active = activate(&first, &self.params.first_scale, &self.params.first_shift);
        let mut second = Array2::zeros((rows, HIDDEN));
        product(
            1.0,
            first_active.view(),
            self.params.second.t(),
            0.0,
            second.view_mut(),
        );
        let second_inverse = normalise(&mut second, &mut self.running[1]);
        let second_active = activate(
            &second,
            &self.params.second_scale,
            &self.params.second_shift,
        );
        let mut output = Array2::from_shape_fn((rows, WIDTH), |(_, at)| self.params.last_bias[at]);
        product(
            1.0,
            second_active.view(),
            self.params.last.t(),
            1.0,
            output.view_mut(),
        );
        let lengths = unit(&mut output);
        Pass {
            input,
            normed: [first, second],
            inverse: [first_inverse, second_inverse],
            active: [first_active, second_active],
            output,
            lengths,
        }
    }

    pub fn backward(&self, pass: &Pass, mut gradient: Array2<f32>) -> Params {
        let rows = gradient.nrows();
        let mut grads = self.params.zeros_like();
        let output = pass.output.as_slice().expect("tables are contiguous");
        rows_mut(&mut gradient)
            .enumerate()
            .for_each(|(at, values)| {
                let unit = &output[at * WIDTH..(at + 1) * WIDTH];
                let along = values.iter().zip(unit).map(|(g, e)| g * e).sum::<f32>();
                for (value, e) in values.iter_mut().zip(unit) {
                    *value = (*value - e * along) / pass.lengths[at];
                }
            });
        product(
            1.0,
            gradient.t(),
            pass.active[1].view(),
            0.0,
            grads.last.view_mut(),
        );
        grads.last_bias = sums(&gradient, None);
        let mut hidden = Array2::zeros((rows, HIDDEN));
        product(
            1.0,
            gradient.view(),
            self.params.last.view(),
            0.0,
            hidden.view_mut(),
        );
        (grads.second_scale, grads.second_shift) =
            self.block_backward(&mut hidden, pass, 1, &self.params.second_scale);
        product(
            1.0,
            hidden.t(),
            pass.active[0].view(),
            0.0,
            grads.second.view_mut(),
        );
        let mut below = Array2::zeros((rows, HIDDEN));
        product(
            1.0,
            hidden.view(),
            self.params.second.view(),
            0.0,
            below.view_mut(),
        );
        (grads.first_scale, grads.first_shift) =
            self.block_backward(&mut below, pass, 0, &self.params.first_scale);
        product(
            1.0,
            below.t(),
            pass.input.view(),
            0.0,
            grads.first.view_mut(),
        );
        grads
    }

    fn block_backward(
        &self,
        gradient: &mut Array2<f32>,
        pass: &Pass,
        at: usize,
        scale: &Array1<f32>,
    ) -> (Array1<f32>, Array1<f32>) {
        let rows = gradient.nrows();
        let active = pass.active[at].as_slice().expect("tables are contiguous");
        rows_mut(gradient).enumerate().for_each(|(row, values)| {
            let active = &active[row * HIDDEN..(row + 1) * HIDDEN];
            for (value, active) in values.iter_mut().zip(active) {
                if *active <= 0.0 {
                    *value = 0.0;
                }
            }
        });
        let normed = &pass.normed[at];
        let shift = sums(gradient, None);
        let scaled = sums(gradient, Some(normed));
        let inverse = &pass.inverse[at];
        let count = rows as f32;
        let normed = normed.as_slice().expect("tables are contiguous");
        rows_mut(gradient).enumerate().for_each(|(row, values)| {
            let normed = &normed[row * HIDDEN..(row + 1) * HIDDEN];
            for (column, value) in values.iter_mut().enumerate() {
                *value = scale[column] * inverse[column] / count
                    * (count * *value - shift[column] - normed[column] * scaled[column]);
            }
        });
        (scaled, shift)
    }

    pub fn embed(&self, input: ArrayView2<f32>) -> Array2<f32> {
        let first = self.layer(input, &self.params.first, 0);
        let second = self.layer(first.view(), &self.params.second, 1);
        let mut output =
            Array2::from_shape_fn((input.nrows(), WIDTH), |(_, at)| self.params.last_bias[at]);
        product(
            1.0,
            second.view(),
            self.params.last.t(),
            1.0,
            output.view_mut(),
        );
        unit(&mut output);
        output
    }

    fn layer(&self, input: ArrayView2<f32>, weights: &Array2<f32>, at: usize) -> Array2<f32> {
        let mut hidden = Array2::zeros((input.nrows(), HIDDEN));
        product(1.0, input, weights.t(), 0.0, hidden.view_mut());
        let (scale, shift) = match at {
            0 => (&self.params.first_scale, &self.params.first_shift),
            _ => (&self.params.second_scale, &self.params.second_shift),
        };
        let running = &self.running[at];
        let inverse = running
            .variance
            .iter()
            .map(|variance| (1.0 / (f64::from(*variance) + EPSILON).sqrt()) as f32)
            .collect::<Vec<_>>();
        rows_mut(&mut hidden).for_each(|values| {
            for (column, value) in values.iter_mut().enumerate() {
                let normed = (*value - running.mean[column]) * inverse[column];
                *value = (normed * scale[column] + shift[column]).max(0.0);
            }
        });
        hidden
    }
}

fn sums(table: &Array2<f32>, times: Option<&Array2<f32>>) -> Array1<f32> {
    let (rows, columns) = table.dim();
    column_sums(rows, columns, |at, sums| match times {
        None => sums
            .iter_mut()
            .zip(row(table, at))
            .for_each(|(sum, value)| *sum += f64::from(*value)),
        Some(times) => sums
            .iter_mut()
            .zip(row(table, at))
            .zip(row(times, at))
            .for_each(|((sum, value), times)| *sum += f64::from(value * times)),
    })
    .into_iter()
    .map(|sum| sum as f32)
    .collect()
}
