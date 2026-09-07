use rand::{Rng, SeedableRng, rngs::StdRng};

use crate::embedding::{features::ContigFeatures, metrics::AggregateMetric};
use crate::refine::bar::MIN_SPLIT_CONTIGS;
use crate::refine::bin_stats::{Centroid, EXACT_LIMIT, centroid};
use crate::refine::dip;

const MAX_ROUNDS: usize = 10;
const FAMILY_ALPHA: f64 = 0.05;

pub fn eligible(features: &ContigFeatures, indices: &[usize], min_bin_size: usize) -> bool {
    indices.len() >= MIN_SPLIT_CONTIGS && features.bin_size(indices) >= 2 * min_bin_size
}

/// Bonferroni over every bin tested this round, with the bootstrap sized so that one draw
/// resolves the corrected level.
pub fn draws_for(eligible: usize) -> usize {
    (eligible.max(1) as f64 / FAMILY_ALPHA).ceil() as usize
}

pub fn candidate(
    features: &ContigFeatures,
    indices: &[usize],
    min_bin_size: usize,
    eligible: usize,
    seed: u64,
) -> Option<[Vec<usize>; 2]> {
    let metric = AggregateMetric::new(features.n_samples() * 2, features.distance_settings())
        .with_bands(features.bands());
    let rows = features.rows(indices);
    let floors = indices
        .iter()
        .map(|index| features.variance_floor(*index))
        .collect::<Vec<_>>();
    let to = |centre: &Centroid| {
        rows.iter()
            .zip(&floors)
            .map(|(row, floor)| metric.distance(row, &centre.row, *floor, centre.floor))
            .collect::<Vec<_>>()
    };

    let from_whole = to(&centroid(features, indices));
    let near = extreme(&from_whole, |a, b| a < b);
    let far = extreme(&from_whole, |a, b| a > b);
    if near == far {
        return None;
    }
    let mut side = rows
        .iter()
        .zip(&floors)
        .map(|(row, floor)| {
            metric.distance(row, &rows[far], *floor, floors[far])
                < metric.distance(row, &rows[near], *floor, floors[near])
        })
        .collect::<Vec<_>>();

    let mut pieces = members(indices, &side)?;
    let mut first = centroid(features, &pieces[0]);
    let mut second = centroid(features, &pieces[1]);
    let mut to_first = to(&first);
    let mut to_second = to(&second);
    for _ in 0..MAX_ROUNDS {
        let next = to_first
            .iter()
            .zip(&to_second)
            .map(|(first, second)| second < first)
            .collect::<Vec<_>>();
        if next == side {
            break;
        }
        side = next;
        pieces = members(indices, &side)?;
        first = centroid(features, &pieces[0]);
        second = centroid(features, &pieces[1]);
        to_first = to(&first);
        to_second = to(&second);
    }
    if pieces
        .iter()
        .any(|piece| features.bin_size(piece) < min_bin_size)
    {
        return None;
    }

    if !bimodal(
        &metric, &to_first, &to_second, &first, &second, eligible, seed,
    ) {
        return None;
    }
    Some(pieces)
}

/// Whether the pieces a split proposes are two modes of the bin rather than two halves of one
/// cloud. Any cut makes tighter pieces, so tightness alone accepts a cut through a genome.
pub fn separates(
    features: &ContigFeatures,
    pieces: &[Vec<usize>],
    eligible: usize,
    seed: u64,
) -> bool {
    let mut order = pieces.iter().collect::<Vec<_>>();
    order.sort_unstable_by_key(|piece| std::cmp::Reverse(features.bin_size(piece)));
    let (Some(largest), Some(next)) = (order.first(), order.get(1)) else {
        return false;
    };
    let indices = pieces.concat();
    let metric = AggregateMetric::new(features.n_samples() * 2, features.distance_settings())
        .with_bands(features.bands());
    let rows = features.rows(&indices);
    let floors = indices
        .iter()
        .map(|index| features.variance_floor(*index))
        .collect::<Vec<_>>();
    let to = |centre: &Centroid| {
        rows.iter()
            .zip(&floors)
            .map(|(row, floor)| metric.distance(row, &centre.row, *floor, centre.floor))
            .collect::<Vec<_>>()
    };
    let first = centroid(features, largest);
    let second = centroid(features, next);
    bimodal(
        &metric,
        &to(&first),
        &to(&second),
        &first,
        &second,
        eligible,
        seed,
    )
}

fn bimodal(
    metric: &AggregateMetric,
    to_first: &[f64],
    to_second: &[f64],
    first: &Centroid,
    second: &Centroid,
    eligible: usize,
    seed: u64,
) -> bool {
    // The difference of two distances saturates at the centroid gap past either centroid,
    // so it piles any cloud up at both ends. The coordinate along the axis does not.
    let gap = metric.distance(&first.row, &second.row, first.floor, second.floor);
    if gap <= 0.0 {
        return false;
    }
    let projection = to_first
        .iter()
        .zip(to_second)
        .map(|(first, second)| (first * first - second * second) / (2.0 * gap))
        .collect::<Vec<_>>();
    // Bases are not observations: a genome in five long contigs is five draws from the
    // mixture, and weighting by length only shrinks the sample the null is drawn from.
    let weights = vec![1.0; projection.len()];
    let (projection, weights) = tested(projection, weights, seed);
    dip::exceeds_null(&projection, &weights, draws_for(eligible), seed)
}

fn extreme(values: &[f64], better: impl Fn(f64, f64) -> bool) -> usize {
    let mut best = 0;
    for (position, value) in values.iter().enumerate() {
        if better(*value, values[best]) {
            best = position;
        }
    }
    best
}

fn members(indices: &[usize], side: &[bool]) -> Option<[Vec<usize>; 2]> {
    let mut pieces = [Vec::new(), Vec::new()];
    for (index, on_second) in indices.iter().zip(side) {
        pieces[usize::from(*on_second)].push(*index);
    }
    if pieces.iter().any(|piece| piece.is_empty()) {
        return None;
    }
    Some(pieces)
}

fn tested(projection: Vec<f64>, weights: Vec<f64>, seed: u64) -> (Vec<f64>, Vec<f64>) {
    let n = projection.len();
    if n <= EXACT_LIMIT {
        return (projection, weights);
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let mut positions = (0..n).collect::<Vec<_>>();
    for position in 0..EXACT_LIMIT {
        positions.swap(position, rng.random_range(position..n));
    }
    positions.truncate(EXACT_LIMIT);
    (
        positions.iter().map(|p| projection[*p]).collect(),
        positions.iter().map(|p| weights[*p]).collect(),
    )
}
