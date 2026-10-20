use std::collections::HashMap;

use ndarray::{Array2, ArrayView2};
use rayon::prelude::*;

use crate::learned::gemm::{fresh, product};

// Each positive pair scored against the negatives alone, as SimCLR's multi-view cross entropy
// (Chen et al. 2020) runs it, so a row's other views never compete with each other.
pub fn per_positive(
    rows: ArrayView2<f32>,
    labels: &[u32],
    temperature: f32,
) -> (f64, Array2<f32>, Vec<(f64, usize)>) {
    let n = rows.nrows();
    let members = members(labels);
    let mut weights = similarity(rows, temperature);
    let losses = weights
        .as_slice_mut()
        .expect("a fresh table is contiguous")
        .par_chunks_mut(n)
        .enumerate()
        .map(|(at, row)| {
            let mates = &members[&labels[at]];
            row[at] = f32::NEG_INFINITY;
            let peak = row.iter().copied().fold(f32::NEG_INFINITY, f32::max);
            let own = mates
                .iter()
                .filter(|mate| **mate != at)
                .map(|mate| (*mate, f64::from(row[*mate])))
                .collect::<Vec<_>>();
            row.iter_mut().for_each(|value| *value = exp(*value - peak));
            mates.iter().for_each(|mate| row[*mate] = 0.0);
            let log_negatives = f64::from(peak) + lane_sum(row).ln();
            let (mut loss, mut spread) = (0.0, 0.0);
            let shares = own
                .iter()
                .map(|(mate, own)| {
                    let high = own.max(log_negatives);
                    let log_total = high + ((own - high).exp() + (log_negatives - high).exp()).ln();
                    loss += log_total - own;
                    spread += (f64::from(peak) - log_total).exp();
                    (*mate, ((own - log_total).exp() - 1.0) as f32)
                })
                .collect::<Vec<_>>();
            let spread = spread as f32;
            row.iter_mut().for_each(|value| *value *= spread);
            shares.iter().for_each(|(mate, share)| row[*mate] = *share);
            (loss, own.len())
        })
        .collect::<Vec<_>>();
    let pairs = losses.iter().map(|(_, pairs)| pairs).sum::<usize>().max(1);
    let loss = losses.iter().map(|(loss, _)| loss).sum::<f64>() / pairs as f64;
    let scale = 1.0 / (temperature * pairs as f32);
    (loss, both_sides(scale, weights, rows), losses)
}

fn members(labels: &[u32]) -> HashMap<u32, Vec<usize>> {
    let mut members = HashMap::<u32, Vec<usize>>::new();
    for (at, label) in labels.iter().enumerate() {
        members.entry(*label).or_default().push(at);
    }
    members
}

// Accelerate runs rows against their own transpose at a tenth of its rate unless the transpose
// is laid out first.
fn similarity(rows: ArrayView2<f32>, temperature: f32) -> Array2<f32> {
    let turned = rows.t().as_standard_layout().into_owned();
    fresh(1.0 / temperature, rows, turned.view())
}

fn both_sides(scale: f32, weights: Array2<f32>, rows: ArrayView2<f32>) -> Array2<f32> {
    let mut gradient = fresh(scale, weights.view(), rows);
    product(scale, weights.t(), rows, 1.0, gradient.view_mut());
    gradient
}

const LANES: usize = 16;

fn lane_sum(values: &[f32]) -> f64 {
    let mut lanes = [0.0f32; LANES];
    let mut chunks = values.chunks_exact(LANES);
    for chunk in &mut chunks {
        lanes
            .iter_mut()
            .zip(chunk)
            .for_each(|(lane, value)| *lane += value);
    }
    let tail = chunks
        .remainder()
        .iter()
        .map(|value| f64::from(*value))
        .sum::<f64>();
    lanes.iter().map(|lane| f64::from(*lane)).sum::<f64>() + tail
}

// A polynomial exp (Moshier 1992, Cephes expf) that the compiler vectorises, where the library
// call runs one value at a time.
fn exp(x: f32) -> f32 {
    let x = x.clamp(-87.0, 88.0);
    let k = (x * std::f32::consts::LOG2_E).round();
    let r = x - k * 0.693_359_4 + k * 2.121_944_4e-4;
    let p = 1.987_569_1e-4_f32;
    let p = p * r + 1.398_199_9e-3;
    let p = p * r + 8.333_452e-3;
    let p = p * r + 4.166_579_6e-2;
    let p = p * r + 1.666_666_5e-1;
    let p = p * r + 0.5;
    let y = p * r * r + r + 1.0;
    y * f32::from_bits(((k as i32 + 127) << 23) as u32)
}
