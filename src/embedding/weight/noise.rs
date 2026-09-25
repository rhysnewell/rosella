use rand::{Rng, rngs::StdRng};
use rayon::prelude::*;

use crate::embedding::metrics::rho;

const RANDOM_PAIRS: usize = 20_000;
pub(crate) const DENSITY_FLOOR: f64 = 1e-12;
const LOG_DEPTH_FLOOR: f64 = 1e-9;
const GRID: usize = 2_048;
const KERNEL_REACH: f64 = 6.0;
const SPREADS: usize = 60;
const SHARES: usize = 49;

// Two halves of one contig share one depth exactly, so each half's partner takes its depth moved
// by the difference to a composition neighbour, weighted by the chance they share a genome.
pub fn partners(
    coverage: &[&[f64]],
    whole: &[&[f64]],
    k: usize,
    rng: &mut StdRng,
) -> Option<Vec<Partners>> {
    let n = coverage.len();
    let samples = coverage[0].len() / 2;
    let depth = |contig: usize, sample: usize| coverage[contig][2 * sample];
    let on = |contig: usize, sample: usize| depth(contig, sample) > 0.0;
    let log_depth = |contig: usize, sample: usize| depth(contig, sample).max(LOG_DEPTH_FLOOR).ln();
    let neighbours = nearest(whole, k);

    let mut log_same = vec![0.0; n * k];
    let mut log_other = vec![0.0; n * k];
    let mut shares = Vec::new();
    for sample in 0..samples {
        let pairs = (0..n * k)
            .filter(|slot| on(slot / k, sample) && on(neighbours[*slot], sample))
            .map(|slot| {
                (
                    slot,
                    log_depth(slot / k, sample) - log_depth(neighbours[slot], sample),
                )
            })
            .collect::<Vec<_>>();
        let mut background = Vec::with_capacity(2 * RANDOM_PAIRS);
        for _ in 0..RANDOM_PAIRS {
            let (a, b) = (rng.random_range(0..n), rng.random_range(0..n));
            if a != b && on(a, sample) && on(b, sample) {
                let difference = log_depth(a, sample) - log_depth(b, sample);
                background.extend([difference, -difference]);
            }
        }
        let Some(density) = Density::fit(&background) else {
            continue;
        };
        if pairs.is_empty() {
            continue;
        }
        let differences = pairs.iter().map(|(_, x)| *x).collect::<Vec<_>>();
        let other = differences
            .iter()
            .map(|x| density.at(*x).max(DENSITY_FLOOR))
            .collect::<Vec<_>>();
        let (spread, share) = fit(&differences, &other);
        shares.push(share);
        for ((slot, x), g) in pairs.iter().zip(&other) {
            log_same[*slot] += log_same_genome(*x, spread);
            log_other[*slot] += g.ln();
        }
    }
    if shares.is_empty() {
        return None;
    }
    let prior = shares.iter().sum::<f64>() / shares.len() as f64;
    let odds = (1.0 - prior).ln() - prior.ln();

    let partners = (0..n)
        .map(|contig| {
            let weights = (0..k)
                .map(|at| {
                    let slot = contig * k + at;
                    1.0 / (1.0
                        + (log_other[slot] + odds - log_same[slot])
                            .clamp(-50.0, 50.0)
                            .exp())
                })
                .collect::<Vec<_>>();
            let total = weights.iter().sum::<f64>();
            let weights = match total > 0.0 {
                true => weights.iter().map(|weight| weight / total).collect(),
                false => vec![1.0 / k as f64; k],
            };
            let rows = (0..k)
                .map(|at| moved(coverage[contig], coverage[neighbours[contig * k + at]]))
                .collect();
            Partners { rows, weights }
        })
        .collect();
    Some(partners)
}

pub struct Partners {
    pub rows: Vec<Vec<f64>>,
    pub weights: Vec<f64>,
}

fn moved(row: &[f64], neighbour: &[f64]) -> Vec<f64> {
    let mut partner = row.to_vec();
    for (sample, (own, other)) in row
        .chunks_exact(2)
        .zip(neighbour.chunks_exact(2))
        .enumerate()
    {
        let (mean, variance) = (own[0], own[1]);
        if mean <= 0.0 {
            continue;
        }
        let shift = match other[0] > 0.0 {
            true => mean.max(LOG_DEPTH_FLOOR).ln() - other[0].max(LOG_DEPTH_FLOOR).ln(),
            false => 0.0,
        };
        let depth = mean * shift.exp();
        partner[2 * sample] = depth;
        partner[2 * sample + 1] = match (variance - mean).abs() <= 1e-8 + 1e-5 * mean.abs() {
            true => depth,
            false => variance * (depth / mean).powi(2),
        };
    }
    partner
}

pub(crate) fn nearest(whole: &[&[f64]], k: usize) -> Vec<usize> {
    (0..whole.len())
        .into_par_iter()
        .flat_map_iter(|contig| {
            let mut others = (0..whole.len())
                .filter(|other| *other != contig)
                .map(|other| (rho(whole[contig], whole[other]), other))
                .collect::<Vec<_>>();
            others.select_nth_unstable_by(k - 1, |a, b| a.0.total_cmp(&b.0));
            others.truncate(k);
            others.into_iter().map(|(_, other)| other)
        })
        .collect()
}

fn log_same_genome(x: f64, spread: f64) -> f64 {
    let variance = 2.0 * spread * spread;
    -x * x / (2.0 * variance) - 0.5 * (2.0 * std::f64::consts::PI * variance).ln()
}

fn fit(differences: &[f64], other: &[f64]) -> (f64, f64) {
    let spreads = (0..SPREADS)
        .map(|at| {
            let (low, high) = (0.01f64.ln(), 3.0f64.ln());
            (low + (high - low) * at as f64 / (SPREADS - 1) as f64).exp()
        })
        .collect::<Vec<_>>();
    let best = |a: (f64, f64, f64), b: (f64, f64, f64)| match b.0 > a.0 {
        true => b,
        false => a,
    };
    let worst = (f64::NEG_INFINITY, 0.0, 0.0);
    let (_, spread, share) = spreads
        .par_iter()
        .map(|spread| {
            let same = differences
                .iter()
                .map(|x| log_same_genome(*x, *spread).exp())
                .collect::<Vec<_>>();
            (0..SHARES)
                .map(|at| {
                    let share = 0.02 + 0.96 * at as f64 / (SHARES - 1) as f64;
                    let likelihood = same
                        .iter()
                        .zip(other)
                        .map(|(f, g)| (share * f + (1.0 - share) * g).ln())
                        .sum::<f64>();
                    (likelihood, *spread, share)
                })
                .fold(worst, best)
        })
        .reduce(|| worst, best);
    (spread, share)
}

// Evaluated once on a grid because every neighbour pair queries it and the background holds
// tens of thousands of points.
pub(crate) struct Density {
    start: f64,
    step: f64,
    values: Vec<f64>,
}

impl Density {
    pub(crate) fn fit(points: &[f64]) -> Option<Self> {
        if points.len() < 2 {
            return None;
        }
        let n = points.len() as f64;
        let mean = points.iter().sum::<f64>() / n;
        let variance = points.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);
        let bandwidth = variance.sqrt() * n.powf(-0.2);
        if !(bandwidth > 0.0) {
            return None;
        }
        let low = points.iter().copied().fold(f64::INFINITY, f64::min) - KERNEL_REACH * bandwidth;
        let high =
            points.iter().copied().fold(f64::NEG_INFINITY, f64::max) + KERNEL_REACH * bandwidth;
        let step = (high - low) / (GRID - 1) as f64;

        let mut mass = vec![0.0; GRID];
        for x in points {
            let at = (x - low) / step;
            let left = (at.floor() as usize).min(GRID - 2);
            let right = at - left as f64;
            mass[left] += 1.0 - right;
            mass[left + 1] += right;
        }
        let reach = (KERNEL_REACH * bandwidth / step).ceil() as usize;
        let kernel = (0..=reach)
            .map(|offset| (-0.5 * (offset as f64 * step / bandwidth).powi(2)).exp())
            .collect::<Vec<_>>();
        let scale = 1.0 / (n * bandwidth * (2.0 * std::f64::consts::PI).sqrt());
        let values = (0..GRID)
            .into_par_iter()
            .map(|at| {
                let first = at.saturating_sub(reach);
                let last = (at + reach).min(GRID - 1);
                (first..=last)
                    .map(|other| mass[other] * kernel[at.abs_diff(other)])
                    .sum::<f64>()
                    * scale
            })
            .collect();
        Some(Self {
            start: low,
            step,
            values,
        })
    }

    pub(crate) fn at(&self, x: f64) -> f64 {
        let at = (x - self.start) / self.step;
        if !(at >= 0.0) || at >= (GRID - 1) as f64 {
            return 0.0;
        }
        let left = at.floor() as usize;
        let right = at - left as f64;
        self.values[left] * (1.0 - right) + self.values[left + 1] * right
    }
}
