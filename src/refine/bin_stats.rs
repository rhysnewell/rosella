use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::embedding::{
    features::ContigFeatures,
    metrics::{euclidean, metabat, rho},
};

pub const METABAT: usize = 0;
pub const RHO: usize = 1;
pub const EUCLIDEAN: usize = 2;
pub const AGGREGATE: usize = 3;

/// Bins above this contribute to the cross-bin thresholds. flight's
/// `min_bin_size_for_averages`.
pub const LARGE_BIN: usize = 1_000_000;

/// All pairs up to here. Past it every contig is scored against one shared sample instead,
/// which keeps the per-contig figures usable where sampling pairs would leave most contigs
/// with no estimate at all.
const EXACT_LIMIT: usize = 2_000;
const REFERENCE_SAMPLE: usize = 1_000;

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

    let weight = features.weight();
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
                let md = metabat(coverage, features.coverage_row(other_index));
                let proportionality = rho(tnf, features.tnf_row(other_index));
                totals[METABAT] += md;
                totals[RHO] += proportionality;
                totals[EUCLIDEAN] += euclidean(tnf, features.tnf_row(other_index));
                totals[AGGREGATE] += (md.powf(weight) * proportionality.powf(1.0 - weight)).sqrt();
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

/// The cross-bin levels a single bin is judged against. flight's `average_bin_stats`.
#[derive(Debug, Clone, Copy, Default)]
pub struct Thresholds {
    pub mean: [f64; 4],
}

impl Thresholds {
    pub fn from_bins<'a>(bins: impl Iterator<Item = (usize, &'a BinStats)>) -> Self {
        let mut totals = [0.0f64; 4];
        let mut counted = 0usize;
        for (bin_size, stats) in bins {
            if bin_size <= LARGE_BIN {
                continue;
            }
            for column in 0..4 {
                totals[column] += stats.mean[column];
            }
            counted += 1;
        }

        if counted == 0 {
            return Self::default();
        }
        Self {
            mean: totals.map(|total| total / counted as f64),
        }
    }
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
    if n <= EXACT_LIMIT {
        return References::All;
    }

    let mut rng = StdRng::seed_from_u64(seed);
    let mut positions = (0..n).collect::<Vec<_>>();
    for position in 0..REFERENCE_SAMPLE {
        positions.swap(position, rng.random_range(position..n));
    }
    positions.truncate(REFERENCE_SAMPLE);
    positions.sort_unstable();
    References::Sampled(positions)
}
