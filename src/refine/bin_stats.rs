use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::embedding::{
    features::ContigFeatures,
    metrics::{euclidean, metabat_with, weight_for},
};
use crate::refine::bar::MIN_SPLIT_CONTIGS;

pub const METABAT: usize = 0;
pub const RHO: usize = 1;
pub const EUCLIDEAN: usize = 2;
pub const AGGREGATE: usize = 3;

/// Bins above this contribute to the cross-bin thresholds. flight's
/// `min_bin_size_for_averages`.
pub const LARGE_BIN: usize = 1_000_000;

pub const SPLIT_LEVEL_NAMES: [&str; 2] = ["flight", "derived"];

/// `Derived` reads a quantile of the run's own spread, so a level follows the assembly rather
/// than a constant carried over from a distance scale this build no longer uses.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum LevelSource {
    #[default]
    Flight,
    Derived,
}

impl LevelSource {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "flight" => Some(Self::Flight),
            "derived" => Some(Self::Derived),
            _ => None,
        }
    }
}

/// All pairs up to here. Past it every contig is scored against one shared sample instead,
/// which keeps the per-contig figures usable where sampling pairs would leave most contigs
/// with no estimate at all.
pub const EXACT_LIMIT: usize = 2_000;
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

    let settings = features.distance_settings();
    let aggregation = settings.aggregation;
    let combination = settings.combination;
    let floors = indices
        .iter()
        .map(|index| features.variance_floor(*index))
        .collect::<Vec<_>>();
    let bands = features.bands();
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
                    aggregation,
                    settings.presence_fraction,
                    bands,
                );
                let proportionality = settings.composition.distance(
                    tnf,
                    features.tnf_row(other_index),
                    settings.composition_scale,
                );
                let weight = weight_for(scored, settings.aggregate_weight);
                totals[METABAT] += md;
                totals[RHO] += proportionality;
                totals[EUCLIDEAN] += euclidean(tnf, features.tnf_row(other_index));
                totals[AGGREGATE] += combination.combine(md, proportionality, weight);
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
        floor += features.variance_floor(*index) * weight;
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

/// The cross-bin levels a single bin is judged against. flight's `average_bin_stats`.
#[derive(Debug, Clone, Copy, Default)]
pub struct Thresholds {
    pub mean: [f64; 4],
    pub source: LevelSource,
}

impl Thresholds {
    pub fn from_bins<'a>(
        bins: impl Iterator<Item = (usize, &'a BinStats)>,
        source: LevelSource,
        quantile: f64,
    ) -> Self {
        let mean = match source {
            LevelSource::Flight => large_bin_means(bins),
            LevelSource::Derived => splittable_quantiles(bins, quantile),
        };
        Self { mean, source }
    }
}

fn large_bin_means<'a>(bins: impl Iterator<Item = (usize, &'a BinStats)>) -> [f64; 4] {
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
        return [0.0; 4];
    }
    totals.map(|total| total / counted as f64)
}

/// Read off the bins that could be split rather than the large ones: a level derived from a
/// population the test never sees describes a different run to the one being judged.
fn splittable_quantiles<'a>(
    bins: impl Iterator<Item = (usize, &'a BinStats)>,
    quantile: f64,
) -> [f64; 4] {
    let mut columns: [Vec<f64>; 4] = Default::default();
    for (_, stats) in bins {
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
