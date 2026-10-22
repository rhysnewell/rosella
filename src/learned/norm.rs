use ndarray::{Array1, Array2};
use rayon::prelude::*;

pub(crate) const EPSILON: f64 = 1e-5;
pub(crate) const MOMENTUM: f32 = 0.1;
// Fixed blocks keep every column sum in one order whatever the thread count.
pub(crate) const BLOCK: usize = 256;

pub struct Running {
    mean: Array1<f32>,
    variance: Array1<f32>,
}

impl Running {
    pub fn new(width: usize) -> Self {
        Self {
            mean: Array1::zeros(width),
            variance: Array1::ones(width),
        }
    }

    pub(crate) fn parts(&self) -> [&[f32]; 2] {
        [&self.mean, &self.variance].map(|values| {
            values
                .as_slice()
                .expect("running statistics are contiguous")
        })
    }

    pub(crate) fn parts_mut(&mut self) -> [&mut [f32]; 2] {
        [&mut self.mean, &mut self.variance].map(|values| {
            values
                .as_slice_mut()
                .expect("running statistics are contiguous")
        })
    }
}

pub struct Normed {
    pub normed: Array2<f32>,
    pub inverse: Array1<f32>,
    pub active: Array2<f32>,
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

fn flat(values: &Array1<f32>) -> &[f32] {
    values.as_slice().expect("vectors are contiguous")
}

// Batch norm (Ioffe & Szegedy 2015) keeps PyTorch's unbiased running variance for embedding.
// Dropout (Srivastava et al. 2014) scales up in training so embedding runs it unchanged.
pub fn normalise(
    mut table: Array2<f32>,
    running: &mut Running,
    (scale, shift): (&Array1<f32>, &Array1<f32>),
    drop: Option<(f32, u64)>,
) -> Normed {
    let (rows, columns) = table.dim();
    let mean = column_sums(rows, columns, |at, sums| {
        sums.iter_mut()
            .zip(row(&table, at))
            .for_each(|(sum, value)| *sum += f64::from(*value));
    })
    .into_iter()
    .map(|sum| sum / rows as f64)
    .collect::<Vec<_>>();
    let variance = column_sums(rows, columns, |at, sums| {
        sums.iter_mut()
            .zip(row(&table, at))
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
    let kept = dropout_scale(drop);
    let bar = dropout_bar(drop);
    let mut active = Array2::<f32>::uninit((rows, columns));
    rows_mut(&mut table)
        .zip(
            active
                .as_slice_mut()
                .expect("a fresh table is contiguous")
                .par_chunks_mut(columns),
        )
        .enumerate()
        .for_each(|(row, (values, active))| {
            let terms = values
                .iter_mut()
                .zip(active.iter_mut())
                .zip(&mean)
                .zip(flat(&inverse))
                .zip(flat(scale))
                .zip(flat(shift));
            match bar {
                None => {
                    for (((((value, active), mean), inverse), scale), shift) in terms {
                        *value = (*value - mean) * inverse;
                        active.write((*value * scale + shift).max(0.0));
                    }
                }
                Some((bar, key)) => {
                    let seed = crate::digest::mix(key ^ row as u64);
                    for (column, (((((value, active), mean), inverse), scale), shift)) in
                        terms.enumerate()
                    {
                        *value = (*value - mean) * inverse;
                        let relu = (*value * scale + shift).max(0.0);
                        let dropped = crate::digest::mix(seed ^ column as u64) < bar;
                        active.write(if dropped { 0.0 } else { relu * kept });
                    }
                }
            }
        });
    Normed {
        normed: table,
        inverse,
        active: unsafe { active.assume_init() },
    }
}

pub fn dropout_scale(drop: Option<(f32, u64)>) -> f32 {
    drop.map_or(1.0, |(share, _)| 1.0 / (1.0 - share))
}

pub(crate) fn dropout_bar(drop: Option<(f32, u64)>) -> Option<(u64, u64)> {
    drop.map(|(share, key)| ((f64::from(share) * u64::MAX as f64) as u64, key))
}

// The mask is a select on bits rather than a branch, since half the units sit at zero at random.
pub fn backward(
    gradient: &mut Array2<f32>,
    step: &Normed,
    scale: &Array1<f32>,
    kept: f32,
) -> (Array1<f32>, Array1<f32>) {
    let (rows, hidden) = gradient.dim();
    let active = step.active.as_slice().expect("tables are contiguous");
    rows_mut(gradient)
        .zip(active.par_chunks(hidden))
        .for_each(|(values, active)| {
            for (value, active) in values.iter_mut().zip(active) {
                let open = u32::from(*active > 0.0).wrapping_neg();
                *value = f32::from_bits((*value * kept).to_bits() & open);
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
    rows_mut(gradient)
        .zip(normed.par_chunks(hidden))
        .for_each(|(values, normed)| {
            for ((((value, normed), factor), shift), scaled) in values
                .iter_mut()
                .zip(normed)
                .zip(&factor)
                .zip(flat(&shift))
                .zip(flat(&scaled))
            {
                *value = factor * (count * *value - shift - normed * scaled);
            }
        });
    (scaled, shift)
}

pub fn sums(table: &Array2<f32>, times: Option<&Array2<f32>>) -> Array1<f32> {
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

pub fn infer(
    mut table: Array2<f32>,
    running: &Running,
    (scale, shift): (&Array1<f32>, &Array1<f32>),
) -> Array2<f32> {
    let inverse = running
        .variance
        .iter()
        .map(|variance| (1.0 / (f64::from(*variance) + EPSILON).sqrt()) as f32)
        .collect::<Vec<_>>();
    rows_mut(&mut table).for_each(|values| {
        for ((((value, mean), inverse), scale), shift) in values
            .iter_mut()
            .zip(flat(&running.mean))
            .zip(&inverse)
            .zip(flat(scale))
            .zip(flat(shift))
        {
            *value = ((*value - mean) * inverse * scale + shift).max(0.0);
        }
    });
    table
}

pub fn unit(table: &mut Array2<f32>) -> Vec<f32> {
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

pub fn unit_backward(gradient: &mut Array2<f32>, output: &Array2<f32>, lengths: &[f32]) {
    let width = output.ncols();
    let output = output.as_slice().expect("tables are contiguous");
    rows_mut(gradient).enumerate().for_each(|(at, values)| {
        let unit = &output[at * width..(at + 1) * width];
        let along = values.iter().zip(unit).map(|(g, e)| g * e).sum::<f32>();
        for (value, e) in values.iter_mut().zip(unit) {
            *value = (*value - e * along) / lengths[at];
        }
    });
}
