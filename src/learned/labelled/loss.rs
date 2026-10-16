use std::collections::HashMap;

use ndarray::{Array2, ArrayView2};
use rayon::prelude::*;

use crate::learned::gemm::product;

// Supervised contrastive loss (Khosla et al. 2020) on unit rows, where a row's positives are
// the other rows carrying its label. A row with no positive anchors nothing.
pub fn contrastive(
    rows: ArrayView2<f32>,
    labels: &[u32],
    temperature: f32,
) -> (f64, Array2<f32>, Vec<(f64, usize)>) {
    let n = rows.nrows();
    let mut members = HashMap::<u32, Vec<usize>>::new();
    for (at, label) in labels.iter().enumerate() {
        members.entry(*label).or_default().push(at);
    }
    let mut weights = Array2::<f32>::zeros((n, n));
    product(1.0 / temperature, rows, rows.t(), 0.0, weights.view_mut());
    let losses = weights
        .as_slice_mut()
        .expect("a fresh table is contiguous")
        .par_chunks_mut(n)
        .enumerate()
        .map(|(at, row)| {
            let mates = &members[&labels[at]];
            let positives = (mates.len() - 1) as f32;
            if positives == 0.0 {
                row.fill(0.0);
                return None;
            }
            row[at] = f32::NEG_INFINITY;
            let peak = row.iter().copied().fold(f32::NEG_INFINITY, f32::max);
            let total = row
                .iter()
                .map(|value| f64::from((value - peak).exp()))
                .sum::<f64>();
            let log_total = f64::from(peak) + total.ln();
            let mated = mates
                .iter()
                .filter(|mate| **mate != at)
                .map(|mate| f64::from(row[*mate]))
                .sum::<f64>();
            let shift = log_total as f32;
            row.iter_mut()
                .for_each(|value| *value = (*value - shift).exp());
            for mate in mates.iter().filter(|mate| **mate != at) {
                row[*mate] -= 1.0 / positives;
            }
            Some(log_total - mated / f64::from(positives))
        })
        .collect::<Vec<_>>();
    let anchors = losses.iter().flatten().count().max(1);
    let loss = losses.iter().flatten().sum::<f64>() / anchors as f64;
    let scale = 1.0 / (temperature * anchors as f32);
    let mut gradient = Array2::zeros(rows.dim());
    product(scale, weights.view(), rows, 0.0, gradient.view_mut());
    product(scale, weights.t(), rows, 1.0, gradient.view_mut());
    let rows = losses
        .iter()
        .map(|loss| loss.map_or((0.0, 0), |loss| (loss, 1)))
        .collect();
    (loss, gradient, rows)
}
