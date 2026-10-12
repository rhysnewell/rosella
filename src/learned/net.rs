use ndarray::{Array1, Array2, ArrayView2, s};
use rand::{Rng, rngs::StdRng};
use rayon::prelude::*;

use crate::learned::gemm::{fresh, later, later_turned, product};

pub const WIDTH: usize = 128;
// Each layer halves the one before. Against a square 2048 second layer this read higher on every
// tier of mouse gut and trained faster.
const LAYERS: [usize; 4] = [2048, 1024, 512, 256];
const EPSILON: f64 = 1e-5;
const MOMENTUM: f32 = 0.1;
// Fixed blocks keep every column sum in one order whatever the thread count.
const BLOCK: usize = 256;

#[derive(Clone)]
struct Layer {
    weights: Array2<f32>,
    scale: Array1<f32>,
    shift: Array1<f32>,
}

#[derive(Clone)]
pub struct Params {
    layers: Vec<Layer>,
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
        let mut from = input;
        let mut layers = Vec::with_capacity(LAYERS.len());
        for width in LAYERS {
            layers.push(Layer {
                weights: uniform(width, from),
                scale: Array1::ones(width),
                shift: Array1::zeros(width),
            });
            from = width;
        }
        let last = uniform(WIDTH, from);
        let last_bias = uniform(1, from).slice(s![0, ..WIDTH]).to_owned();
        Self {
            layers,
            last,
            last_bias,
        }
    }

    pub fn zeros_like(&self) -> Self {
        Self {
            layers: self
                .layers
                .iter()
                .map(|layer| Layer {
                    weights: Array2::zeros(layer.weights.dim()),
                    scale: Array1::zeros(layer.scale.len()),
                    shift: Array1::zeros(layer.shift.len()),
                })
                .collect(),
            last: Array2::zeros(self.last.dim()),
            last_bias: Array1::zeros(WIDTH),
        }
    }

    pub fn slices(&self) -> Vec<&[f32]> {
        let mut out = Vec::new();
        for layer in &self.layers {
            out.extend([flat(&layer.weights), flat(&layer.scale), flat(&layer.shift)]);
        }
        out.extend([flat(&self.last), flat(&self.last_bias)]);
        out
    }

    pub fn slices_mut(&mut self) -> Vec<&mut [f32]> {
        let mut out = Vec::new();
        for layer in &mut self.layers {
            out.extend([
                flat_mut(&mut layer.weights),
                flat_mut(&mut layer.scale),
                flat_mut(&mut layer.shift),
            ]);
        }
        out.extend([flat_mut(&mut self.last), flat_mut(&mut self.last_bias)]);
        out
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
    running: Vec<Running>,
}

struct Step {
    normed: Array2<f32>,
    inverse: Array1<f32>,
    active: Array2<f32>,
}

pub struct Pass {
    input: Array2<f32>,
    steps: Vec<Step>,
    pub output: Array2<f32>,
    lengths: Vec<f32>,
    kept: f32,
}

impl Pass {
    fn below(&self, at: usize) -> ArrayView2<'_, f32> {
        match at {
            0 => self.input.view(),
            _ => self.steps[at - 1].active.view(),
        }
    }
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
// 2015), the running variance unbiased as PyTorch keeps it. The same pass writes the activation,
// with dropout (Srivastava et al. 2014) scaled up in training so embedding runs it unchanged.
fn normalise(
    table: &mut Array2<f32>,
    running: &mut Running,
    (scale, shift): (&Array1<f32>, &Array1<f32>),
    drop: Option<(f32, u64)>,
) -> (Array1<f32>, Array2<f32>, f32) {
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
    let kept = drop.map_or(1.0, |(share, _)| 1.0 / (1.0 - share));
    let bar = drop.map(|(share, key)| ((f64::from(share) * u64::MAX as f64) as u64, key));
    let mut active = Array2::<f32>::uninit((rows, columns));
    rows_mut(table)
        .zip(
            active
                .as_slice_mut()
                .expect("a fresh table is contiguous")
                .par_chunks_mut(columns),
        )
        .enumerate()
        .for_each(|(row, (values, active))| {
            let seed = bar.map(|(bar, key)| (bar, crate::digest::mix(key ^ row as u64)));
            let terms = values
                .iter_mut()
                .zip(active.iter_mut())
                .zip(&mean)
                .zip(&inverse)
                .zip(scale)
                .zip(shift)
                .enumerate();
            for (column, (((((value, active), mean), inverse), scale), shift)) in terms {
                *value = (*value - mean) * inverse;
                let relu = (*value * scale + shift).max(0.0);
                active.write(match seed {
                    Some((bar, seed)) if crate::digest::mix(seed ^ column as u64) < bar => 0.0,
                    Some(_) => relu * kept,
                    None => relu,
                });
            }
        });
    (inverse, unsafe { active.assume_init() }, kept)
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
        let running = LAYERS
            .iter()
            .map(|&width| Running {
                mean: Array1::zeros(width),
                variance: Array1::ones(width),
            })
            .collect();
        Self {
            params: Params::new(input, rng),
            running,
        }
    }

    pub fn forward(&mut self, input: Array2<f32>) -> Pass {
        self.forward_dropping(input, None)
    }

    pub fn forward_dropping(&mut self, input: Array2<f32>, drop: Option<(f32, u64)>) -> Pass {
        let rows = input.nrows();
        let mut steps: Vec<Step> = Vec::with_capacity(LAYERS.len());
        let mut kept = 1.0;
        for (at, (layer, running)) in self.params.layers.iter().zip(&mut self.running).enumerate() {
            let below = steps.last().map_or(input.view(), |step| step.active.view());
            let mut normed = fresh(1.0, below, layer.weights.t());
            let drop = drop.filter(|_| at == 0);
            let (inverse, active, share) =
                normalise(&mut normed, running, (&layer.scale, &layer.shift), drop);
            if at == 0 {
                kept = share;
            }
            steps.push(Step {
                normed,
                inverse,
                active,
            });
        }
        let top = steps.last().map_or(input.view(), |step| step.active.view());
        let mut output = Array2::from_shape_fn((rows, WIDTH), |(_, at)| self.params.last_bias[at]);
        product(1.0, top, self.params.last.t(), 1.0, output.view_mut());
        let lengths = unit(&mut output);
        Pass {
            input,
            steps,
            output,
            lengths,
            kept,
        }
    }

    pub fn backward(&self, pass: &Pass, gradient: Array2<f32>) -> Params {
        self.backward_parts(pass, gradient, false).0
    }

    pub fn backward_input(&self, pass: &Pass, gradient: Array2<f32>) -> (Params, Array2<f32>) {
        let (grads, input) = self.backward_parts(pass, gradient, true);
        (grads, input.expect("asked for the input gradient"))
    }

    fn backward_parts(
        &self,
        pass: &Pass,
        mut gradient: Array2<f32>,
        want_input: bool,
    ) -> (Params, Option<Array2<f32>>) {
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
        let layers = self.params.layers.len();
        let last = later(1.0, gradient.t(), pass.below(layers));
        let last_bias = sums(&gradient, None);
        let mut below = fresh(1.0, gradient.view(), self.params.last.view());
        let mut grads = Vec::with_capacity(layers);
        let mut input = None;
        for at in (0..layers).rev() {
            let (layer, step) = (&self.params.layers[at], &pass.steps[at]);
            if at == 0 && pass.kept != 1.0 {
                rows_mut(&mut below)
                    .for_each(|values| values.iter_mut().for_each(|v| *v *= pass.kept));
            }
            let (scale, shift) = block_backward(&mut below, step, &layer.scale);
            let next =
                (at > 0 || want_input).then(|| fresh(1.0, below.view(), layer.weights.view()));
            grads.push((later_turned(1.0, below, pass.below(at)), scale, shift));
            match next {
                Some(next) if at > 0 => below = next,
                next => {
                    input = next;
                    break;
                }
            }
        }
        let layers = grads
            .into_iter()
            .rev()
            .map(|(weights, scale, shift)| Layer {
                weights: weights.wait(),
                scale,
                shift,
            })
            .collect();
        let grads = Params {
            layers,
            last: last.wait(),
            last_bias,
        };
        (grads, input)
    }

    pub fn embed(&self, input: ArrayView2<f32>) -> Array2<f32> {
        let mut hidden: Option<Array2<f32>> = None;
        for (layer, running) in self.params.layers.iter().zip(&self.running) {
            let below = hidden.as_ref().map_or(input, |table| table.view());
            let mut raw = fresh(1.0, below, layer.weights.t());
            let inverse = running
                .variance
                .iter()
                .map(|variance| (1.0 / (f64::from(*variance) + EPSILON).sqrt()) as f32)
                .collect::<Vec<_>>();
            rows_mut(&mut raw).for_each(|values| {
                for (column, value) in values.iter_mut().enumerate() {
                    let normed = (*value - running.mean[column]) * inverse[column];
                    *value = (normed * layer.scale[column] + layer.shift[column]).max(0.0);
                }
            });
            hidden = Some(raw);
        }
        let top = hidden.as_ref().map_or(input, |table| table.view());
        let mut output =
            Array2::from_shape_fn((input.nrows(), WIDTH), |(_, at)| self.params.last_bias[at]);
        product(1.0, top, self.params.last.t(), 1.0, output.view_mut());
        unit(&mut output);
        output
    }
}

fn block_backward(
    gradient: &mut Array2<f32>,
    step: &Step,
    scale: &Array1<f32>,
) -> (Array1<f32>, Array1<f32>) {
    let (rows, hidden) = gradient.dim();
    let active = step.active.as_slice().expect("tables are contiguous");
    rows_mut(gradient).enumerate().for_each(|(row, values)| {
        let active = &active[row * hidden..(row + 1) * hidden];
        for (value, active) in values.iter_mut().zip(active) {
            *value = if *active > 0.0 { *value } else { 0.0 };
        }
    });
    let shift = sums(gradient, None);
    let scaled = sums(gradient, Some(&step.normed));
    let count = rows as f32;
    let factor = scale
        .iter()
        .zip(&step.inverse)
        .map(|(scale, inverse)| scale * inverse / count)
        .collect::<Vec<_>>();
    let normed = step.normed.as_slice().expect("tables are contiguous");
    rows_mut(gradient).enumerate().for_each(|(row, values)| {
        let normed = &normed[row * hidden..(row + 1) * hidden];
        for ((((value, normed), factor), shift), scaled) in values
            .iter_mut()
            .zip(normed)
            .zip(&factor)
            .zip(&shift)
            .zip(&scaled)
        {
            *value = factor * (count * *value - shift - normed * scaled);
        }
    });
    (scaled, shift)
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
