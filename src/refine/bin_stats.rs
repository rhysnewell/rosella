use rayon::prelude::*;

use crate::embedding::{
    features::ContigFeatures,
    metrics::{Point, aggregate, euclidean, rho_between},
};
use crate::refine::bar::MIN_SPLIT_CONTIGS;

pub const ABUNDANCE: usize = 0;
pub const RHO: usize = 1;
pub const EUCLIDEAN: usize = 2;
pub const AGGREGATE: usize = 3;

pub struct BinStats {
    pub mean: [f64; 4],
    pub std: [f64; 4],
    pub per_contig: Vec<[f64; 4]>,
}

pub fn bin_stats(features: &ContigFeatures, indices: &[usize], seed: u64) -> Option<BinStats> {
    if indices.len() < 2 {
        return None;
    }

    let points = indices
        .par_iter()
        .map(|index| features.point(*index))
        .collect::<Vec<_>>();
    let pair = Pair {
        features,
        indices,
        points: &points,
        weight: features.distance_settings().aggregate_weight,
    };
    let per_contig = match indices.len() <= crate::tuning::EXACT_LIMIT {
        true => pair.all_pairs(),
        false => pair.against(&references(indices.len(), seed)),
    };

    let mut mean = [0.0f64; 4];
    let mut std = [0.0f64; 4];
    for column in 0..4 {
        let values = per_contig.iter().map(|row| row[column]);
        mean[column] = values.clone().sum::<f64>() / per_contig.len() as f64;
        let variance = values
            .map(|value| (value - mean[column]) * (value - mean[column]))
            .sum::<f64>()
            / per_contig.len() as f64;
        std[column] = variance.sqrt();
    }

    Some(BinStats {
        mean,
        std,
        per_contig,
    })
}

struct Pair<'a> {
    features: &'a ContigFeatures<'a>,
    indices: &'a [usize],
    points: &'a [Point],
    weight: Option<f64>,
}

impl Pair<'_> {
    fn between(&self, a: usize, b: usize) -> [f64; 4] {
        let (left, right) = (&self.points[a], &self.points[b]);
        let proportionality = rho_between(&left.composition, &right.composition);
        let (abundance, _, together) = aggregate(
            &left.abundance,
            &right.abundance,
            proportionality,
            self.weight,
        );
        [
            abundance,
            proportionality,
            euclidean(
                self.features.tnf_row(self.indices[a]),
                self.features.tnf_row(self.indices[b]),
            ),
            together,
        ]
    }

    // Each pair is measured once and added to both ends. Rows run in order, so every contig
    // still sums its partners in index order and the means match a full scan to the bit.
    fn all_pairs(&self) -> Vec<[f64; 4]> {
        let n = self.indices.len();
        let mut totals = vec![[0.0f64; 4]; n];
        let mut row = Vec::with_capacity(n);
        for at in 0..n {
            let (head, tail) = totals.split_at_mut(at + 1);
            tail.par_iter_mut()
                .with_min_len(PAIRS_PER_TASK)
                .enumerate()
                .map(|(offset, total)| {
                    let distances = self.between(at, at + 1 + offset);
                    add(total, &distances);
                    distances
                })
                .collect_into_vec(&mut row);
            for distances in &row {
                add(&mut head[at], distances);
            }
        }
        let divisor = (n - 1) as f64;
        totals
            .into_iter()
            .map(|total| total.map(|sum| sum / divisor))
            .collect()
    }

    fn against(&self, references: &[usize]) -> Vec<[f64; 4]> {
        (0..self.indices.len())
            .into_par_iter()
            .map(|position| {
                let mut totals = [0.0f64; 4];
                let mut counted = 0usize;
                for other in references.iter().filter(|other| **other != position) {
                    add(&mut totals, &self.between(position, *other));
                    counted += 1;
                }
                let divisor = counted.max(1) as f64;
                totals.map(|total| total / divisor)
            })
            .collect()
    }
}

const PAIRS_PER_TASK: usize = 64;

fn add(total: &mut [f64; 4], distances: &[f64; 4]) {
    for (sum, distance) in total.iter_mut().zip(distances) {
        *sum += distance;
    }
}

pub fn centroid(features: &ContigFeatures, indices: &[usize]) -> Point {
    let coverage_columns = features.n_samples() * 2;
    let tnf_columns = features.tnf_row(indices[0]).len();
    let mut row = vec![0.0; coverage_columns + tnf_columns];
    let mut total = 0.0;

    for index in indices {
        let weight = features.length(*index) as f64;
        for (slot, value) in row[..coverage_columns]
            .iter_mut()
            .zip(features.coverage_row(*index))
        {
            *slot += value * weight;
        }
        for (slot, value) in row[coverage_columns..]
            .iter_mut()
            .zip(features.tnf_row(*index))
        {
            *slot += value * weight;
        }
        total += weight;
    }

    for slot in row.iter_mut() {
        *slot /= total;
    }
    let (coverage, composition) = row.split_at(coverage_columns);
    Point::new(coverage, composition, features.distance_settings())
}

// The cross-bin levels a single bin is judged against, read off the bins that could be split.
#[derive(Debug, Clone, Copy, Default)]
pub struct Thresholds {
    pub quantile: [f64; 4],
}

impl Thresholds {
    pub fn from_bins<'a>(bins: impl Iterator<Item = &'a BinStats>) -> Self {
        Self {
            quantile: splittable_quantiles(bins, crate::tuning::SPLIT_LEVEL_QUANTILE),
        }
    }
}

// Read off the bins that could be split rather than the large ones: a level derived from a
// population the test never sees describes a different run to the one being judged.
fn splittable_quantiles<'a>(bins: impl Iterator<Item = &'a BinStats>, quantile: f64) -> [f64; 4] {
    let mut columns: [Vec<f64>; 4] = Default::default();
    for stats in bins {
        if stats.per_contig.len() < MIN_SPLIT_CONTIGS {
            continue;
        }
        for (column, spreads) in columns.iter_mut().enumerate() {
            spreads.push(stats.mean[column]);
        }
    }

    let mut levels = [0.0f64; 4];
    for (level, spreads) in levels.iter_mut().zip(columns.iter_mut()) {
        *level = percentile(spreads, quantile);
    }
    levels
}

fn percentile(values: &mut [f64], quantile: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(f64::total_cmp);
    let position = quantile.clamp(0.0, 1.0) * (values.len() - 1) as f64;
    let below = position.floor() as usize;
    let above = position.ceil() as usize;
    values[below] + (values[above] - values[below]) * (position - below as f64)
}

fn references(n: usize, seed: u64) -> Vec<usize> {
    let mut positions = crate::seeds::sample_positions(n, crate::tuning::REFERENCE_SAMPLE, seed);
    positions.sort_unstable();
    positions
}
