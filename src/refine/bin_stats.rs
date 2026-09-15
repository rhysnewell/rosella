use rayon::prelude::*;

use crate::embedding::{
    features::ContigFeatures,
    metrics::{combine, euclidean, metabat_with, rho, weight_for},
};
use crate::refine::bar::MIN_SPLIT_CONTIGS;

pub const METABAT: usize = 0;
pub const RHO: usize = 1;
pub const EUCLIDEAN: usize = 2;
pub const AGGREGATE: usize = 3;

/// Mean metabat, rho, tetranucleotide euclidean and aggregate distance within a bin, both
/// per contig and across the bin. flight's `metrics.get_averages`.
pub struct BinStats {
    pub mean: [f64; 4],
    pub std: [f64; 4],
    pub per_contig: Vec<[f64; 4]>,
}

pub fn bin_stats(features: &ContigFeatures, indices: &[usize], seed: u64) -> Option<BinStats> {
    if indices.len() < 2 {
        return None;
    }

    let settings = features.distance_settings();
    let floors = features.floors(indices);
    let references = references(indices.len(), seed);

    let per_contig = indices
        .par_iter()
        .enumerate()
        .map(|(position, index)| {
            let coverage = features.coverage_row(*index);
            let tnf = features.tnf_row(*index);

            let mut totals = [0.0f64; 4];
            let mut counted = 0usize;
            for other in references.iter(indices.len()) {
                if other == position {
                    continue;
                }
                let other_index = indices[other];
                let (md, scored) = metabat_with(
                    coverage,
                    features.coverage_row(other_index),
                    floors[position],
                    floors[other],
                    settings.presence_fraction,
                );
                let proportionality = rho(tnf, features.tnf_row(other_index));
                let weight = weight_for(scored, None);
                totals[METABAT] += md;
                totals[RHO] += proportionality;
                totals[EUCLIDEAN] += euclidean(tnf, features.tnf_row(other_index));
                totals[AGGREGATE] += combine(md, proportionality, weight);
                counted += 1;
            }

            let divisor = counted.max(1) as f64;
            [
                totals[METABAT] / divisor,
                totals[RHO] / divisor,
                totals[EUCLIDEAN] / divisor,
                totals[AGGREGATE] / divisor,
            ]
        })
        .collect::<Vec<_>>();

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

pub struct Centroid {
    pub row: Vec<f64>,
    pub floor: f64,
}

pub fn centroid(features: &ContigFeatures, indices: &[usize]) -> Centroid {
    let coverage_columns = features.n_samples() * 2;
    let tnf_columns = features.tnf_row(indices[0]).len();
    let mut row = vec![0.0; coverage_columns + tnf_columns];
    let mut floor = 0.0;
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
        floor += crate::embedding::metrics::MIN_VAR * weight;
        total += weight;
    }

    for slot in row.iter_mut() {
        *slot /= total;
    }
    Centroid {
        row,
        floor: floor / total,
    }
}

/// The cross-bin levels a single bin is judged against, read off the bins that could be split.
#[derive(Debug, Clone, Copy, Default)]
pub struct Thresholds {
    pub mean: [f64; 4],
}

impl Thresholds {
    pub fn from_bins<'a>(bins: impl Iterator<Item = &'a BinStats>, quantile: f64) -> Self {
        Self {
            mean: splittable_quantiles(bins, quantile),
        }
    }
}

/// Read off the bins that could be split rather than the large ones: a level derived from a
/// population the test never sees describes a different run to the one being judged.
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

enum References {
    All,
    Sampled(Vec<usize>),
}

impl References {
    fn iter(&self, n: usize) -> Box<dyn Iterator<Item = usize> + '_> {
        match self {
            References::All => Box::new(0..n),
            References::Sampled(positions) => Box::new(positions.iter().copied()),
        }
    }
}

fn references(n: usize, seed: u64) -> References {
    if n <= crate::tuning::EXACT_LIMIT {
        return References::All;
    }

    let mut positions =
        crate::seeds::sample_positions(n, crate::tuning::REFERENCE_SAMPLE, seed);
    positions.sort_unstable();
    References::Sampled(positions)
}
