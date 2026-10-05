use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::embedding::metrics::likelihood::{CoverageDistance, DepthModel, SampleModel};

use super::noise::{DENSITY_FLOOR, Density, RANDOM_PAIRS, nearest};

pub const QUERIES: usize = 2_000;
pub const POOL: usize = 10_000;
const NEIGHBOURS: usize = 10;
const MIN_PAIRS: usize = 50;
const READ_LENGTHS: usize = 13;
const LOWEST_READ_LENGTH: f64 = 1.5;
const READ_LENGTH_STEP: f64 = 0.25;
const SPREADS: usize = 30;
const LOWEST_SPREAD: f64 = 0.01;
const HIGHEST_SPREAD: f64 = 3.0;
const SHARE_STEPS: usize = 8;
const LOWEST_SHARE: f64 = 0.02;

struct Pair {
    gap: f64,
    counts: f64,
    other: f64,
}

struct Fitted {
    model: SampleModel,
    share: f64,
    log_likelihood: f64,
}

// Composition neighbours stand in for genome mates, so nothing here reads a coverage distance and
// the model can be fitted before the first one is needed. The pool is a random draw of long contigs.
pub fn fit(
    coverage: &[&[f64]],
    lengths: &[f64],
    composition: &[&[f64]],
    form: CoverageDistance,
    seed: u64,
) -> Option<DepthModel> {
    let queries = coverage.len().min(QUERIES);
    if queries <= NEIGHBOURS {
        return None;
    }
    let neighbours = nearest(composition, queries, NEIGHBOURS);
    let samples = coverage[0].len() / 2;
    let fitted = (0..samples)
        .into_par_iter()
        .map(|sample| {
            let depth = |contig: usize| coverage[contig][2 * sample];
            let mut rng = StdRng::seed_from_u64(seed.wrapping_add(sample as u64));
            fit_sample(&depth, lengths, &neighbours, queries, &mut rng)
        })
        .collect::<Vec<_>>();
    let shares = fitted
        .iter()
        .flatten()
        .map(|fitted| fitted.share)
        .collect::<Vec<_>>();
    if shares.is_empty() {
        return None;
    }
    let prior = shares.iter().sum::<f64>() / shares.len() as f64;
    let model = DepthModel {
        form,
        prior: (1.0 - prior).ln() - prior.ln(),
        samples: fitted
            .into_iter()
            .map(|fitted| fitted.map(|fitted| fitted.model))
            .collect(),
    };
    Some(model)
}

fn fit_sample(
    depth: &impl Fn(usize) -> f64,
    lengths: &[f64],
    neighbours: &[usize],
    queries: usize,
    rng: &mut StdRng,
) -> Option<Fitted> {
    let n = lengths.len();
    let present = (0..n).filter(|contig| depth(*contig) > 0.0).count();
    let mut background = Vec::with_capacity(2 * RANDOM_PAIRS);
    for _ in 0..RANDOM_PAIRS {
        let (a, b) = (rng.random_range(0..n), rng.random_range(0..n));
        if a != b && depth(a) > 0.0 && depth(b) > 0.0 {
            let gap = depth(a).ln() - depth(b).ln();
            background.extend([gap, -gap]);
        }
    }
    let density = Density::fit(&background)?;
    let reads = |contig: usize| depth(contig) * lengths[contig];
    let pairs = (0..queries * NEIGHBOURS)
        .map(|slot| (slot / NEIGHBOURS, neighbours[slot]))
        .filter(|(a, b)| depth(*a) > 0.0 && depth(*b) > 0.0)
        .map(|(a, b)| {
            let gap = depth(a).ln() - depth(b).ln();
            Pair {
                gap,
                counts: 1.0 / reads(a) + 1.0 / reads(b),
                other: density.at(gap).max(DENSITY_FLOOR),
            }
        })
        .collect::<Vec<_>>();
    if pairs.len() < MIN_PAIRS {
        return None;
    }
    let mut best: Option<Fitted> = None;
    for read_length in (0..READ_LENGTHS).map(read_length_at) {
        for spread in (0..SPREADS).map(spread_at) {
            let same = pairs
                .iter()
                .map(|pair| normal(pair.gap, read_length * pair.counts + 2.0 * spread))
                .collect::<Vec<_>>();
            let (share, log_likelihood) = mixture(&same, &pairs);
            if best
                .as_ref()
                .is_none_or(|best| log_likelihood > best.log_likelihood)
            {
                best = Some(Fitted {
                    model: SampleModel {
                        read_length,
                        spread,
                        absent: ((n - present + 1) as f64 / (n + 2) as f64).ln(),
                        background: density.logged(),
                    },
                    share,
                    log_likelihood,
                });
            }
        }
    }
    best
}

// The likelihood is concave in the share, so a few Newton steps find it.
fn mixture(same: &[f64], pairs: &[Pair]) -> (f64, f64) {
    let mut share = 0.5;
    for _ in 0..SHARE_STEPS {
        let (slope, curve) =
            same.iter()
                .zip(pairs)
                .fold((0.0, 0.0), |(slope, curve), (f, pair)| {
                    let lean = (f - pair.other) / (share * f + (1.0 - share) * pair.other);
                    (slope + lean, curve + lean * lean)
                });
        if curve <= 0.0 {
            break;
        }
        share = (share + slope / curve).clamp(LOWEST_SHARE, 1.0 - LOWEST_SHARE);
    }
    let log_likelihood = same
        .iter()
        .zip(pairs)
        .map(|(f, pair)| (share * f + (1.0 - share) * pair.other).ln())
        .sum();
    (share, log_likelihood)
}

fn normal(x: f64, variance: f64) -> f64 {
    (-0.5 * x * x / variance).exp() / (std::f64::consts::TAU * variance).sqrt()
}

pub fn read_length_at(step: usize) -> f64 {
    10f64.powf(LOWEST_READ_LENGTH + step as f64 * READ_LENGTH_STEP)
}

fn spread_at(step: usize) -> f64 {
    let (low, high) = (LOWEST_SPREAD.ln(), HIGHEST_SPREAD.ln());
    (low + (high - low) * step as f64 / (SPREADS - 1) as f64)
        .exp()
        .powi(2)
}
