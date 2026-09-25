use ndarray::Array2;
use rand::{Rng, rngs::StdRng};
use rayon::prelude::*;

use super::{LEAST_POINTS, Scatter, candidates};
use crate::embedding::weight::noise::{DENSITY_FLOOR, Density};

const RANDOM_PAIRS: usize = 20_000;
const STEPS: usize = 31;
const SHARE_ROUNDS: usize = 30;
const SHARE_BOUNDS: (f64, f64) = (1e-3, 1.0 - 1e-3);

struct Pair {
    difference: f64,
    reach: f64,
}

// Bins are drawn partly on depth, so their members understate how far a genome's depth strays.
// Composition neighbours are drawn blind to depth, and a mixture against random pairs keeps only
// the share that looks like one genome.
pub fn fit_neighbours(
    table: &Array2<f64>,
    lengths: &[usize],
    pool: &[usize],
    neighbours: &[usize],
    k: usize,
    rng: &mut StdRng,
) -> Vec<Option<Scatter>> {
    (0..table.ncols() / 2)
        .map(|sample| {
            let depth = |at: usize| table[[pool[at], 2 * sample]];
            let reach = |at: usize| 1.0 / (depth(at) * lengths[pool[at]].max(1) as f64);
            let pairs = (0..pool.len() * k)
                .map(|slot| (slot / k, neighbours[slot]))
                .filter(|(a, b)| depth(*a) > 0.0 && depth(*b) > 0.0)
                .map(|(a, b)| Pair {
                    difference: depth(a).ln() - depth(b).ln(),
                    reach: reach(a) + reach(b),
                })
                .collect::<Vec<_>>();
            let mut background = Vec::with_capacity(2 * RANDOM_PAIRS);
            for _ in 0..RANDOM_PAIRS {
                let (a, b) = (
                    rng.random_range(0..pool.len()),
                    rng.random_range(0..pool.len()),
                );
                if a != b && depth(a) > 0.0 && depth(b) > 0.0 {
                    let difference = depth(a).ln() - depth(b).ln();
                    background.extend([difference, -difference]);
                }
            }
            if pairs.len() < LEAST_POINTS {
                return None;
            }
            let density = Density::fit(&background)?;
            let other = pairs
                .iter()
                .map(|pair| density.at(pair.difference).max(DENSITY_FLOOR))
                .collect::<Vec<_>>();
            candidates(STEPS)
                .into_par_iter()
                .map(|model| (likelihood(&model, &pairs, &other), model))
                .max_by(|a, b| a.0.total_cmp(&b.0))
                .map(|(_, model)| model)
        })
        .collect()
}

fn likelihood(model: &Scatter, pairs: &[Pair], other: &[f64]) -> f64 {
    let same = pairs
        .iter()
        .map(|pair| {
            let variance = model.sampling * pair.reach + 2.0 * model.bias;
            (-pair.difference * pair.difference / (2.0 * variance)).exp()
                / (2.0 * std::f64::consts::PI * variance).sqrt()
        })
        .collect::<Vec<_>>();
    let mut share = 0.5;
    for _ in 0..SHARE_ROUNDS {
        share = (same
            .iter()
            .zip(other)
            .map(|(f, g)| share * f / (share * f + (1.0 - share) * g))
            .sum::<f64>()
            / same.len() as f64)
            .clamp(SHARE_BOUNDS.0, SHARE_BOUNDS.1);
    }
    same.iter()
        .zip(other)
        .map(|(f, g)| (share * f + (1.0 - share) * g).ln())
        .sum()
}
