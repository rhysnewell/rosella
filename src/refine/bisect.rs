use crate::embedding::{features::ContigFeatures, metrics::AggregateMetric};
use crate::refine::bar::MIN_SPLIT_CONTIGS;
use crate::refine::bin_stats::{Centroid, centroid};
use crate::refine::dip;

const MAX_ROUNDS: usize = 10;

pub fn eligible(features: &ContigFeatures, indices: &[usize], min_bin_size: usize) -> bool {
    indices.len() >= MIN_SPLIT_CONTIGS
        && features.bin_size(indices) >= crate::tuning::BISECT_SIZE_MULTIPLE * min_bin_size
}

/// Bonferroni over every bin tested this round, with the bootstrap sized so that one draw
/// resolves the corrected level.
fn draws_for(eligible: usize) -> usize {
    (eligible.max(1) as f64 / crate::tuning::FAMILY_ALPHA).ceil() as usize
}

struct Projector {
    metric: AggregateMetric,
    rows: Vec<Vec<f64>>,
    floors: Vec<f64>,
}

impl Projector {
    fn new(features: &ContigFeatures, indices: &[usize]) -> Self {
        Self {
            metric: AggregateMetric::new(features.n_samples() * 2, features.distance_settings()),
            rows: features.rows(indices),
            floors: features.floors(indices),
        }
    }

    fn to(&self, centre: &Centroid) -> Vec<f64> {
        self.rows
            .iter()
            .zip(&self.floors)
            .map(|(row, floor)| self.metric.distance(row, &centre.row, *floor, centre.floor))
            .collect()
    }

    fn nearer(&self, first: usize, second: usize) -> Vec<bool> {
        self.rows
            .iter()
            .zip(&self.floors)
            .map(|(row, floor)| {
                self.metric
                    .distance(row, &self.rows[first], *floor, self.floors[first])
                    < self
                        .metric
                        .distance(row, &self.rows[second], *floor, self.floors[second])
            })
            .collect()
    }
}

/// Whether the bin's own contigs make two clouds, one cloud, or too little to ask. Hartigan and
/// Hartigan (1985) against a uniform null, so the bar is a significance level and not a constant.
enum Shape {
    Untestable(Mute),
    OneCloud,
    TwoClouds([Vec<usize>; 2]),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mute {
    NoAxis,
    OneSided,
    PieceTooSmall,
}

#[derive(Clone, Copy)]
struct Trial {
    min_bin_size: usize,
    eligible: usize,
    seed: u64,
}

fn shape(
    features: &ContigFeatures,
    indices: &[usize],
    min_bin_size: usize,
    eligible: usize,
    seed: u64,
) -> Shape {
    let project = Projector::new(features, indices);
    let from_whole = project.to(&centroid(features, indices));
    let near = extreme(&from_whole, |a, b| a < b);
    let far = extreme(&from_whole, |a, b| a > b);
    if near == far {
        return Shape::Untestable(Mute::NoAxis);
    }
    grow(
        &project,
        features,
        indices,
        (far, near),
        Trial {
            min_bin_size,
            eligible,
            seed,
        },
    )
}

pub fn candidate(
    features: &ContigFeatures,
    indices: &[usize],
    min_bin_size: usize,
    eligible: usize,
    seed: u64,
) -> Option<[Vec<usize>; 2]> {
    match shape(features, indices, min_bin_size, eligible, seed) {
        Shape::TwoClouds(pieces) => Some(pieces),
        _ => None,
    }
}

/// The markers can name the two contigs a fused bin is fused from, which is a better pair to
/// grow from than the two the geometry happens to put furthest apart.
pub fn from_seeds(
    features: &ContigFeatures,
    indices: &[usize],
    seeds: (usize, usize),
    min_bin_size: usize,
    eligible: usize,
    seed: u64,
) -> Option<[Vec<usize>; 2]> {
    if seeds.0 == seeds.1 || seeds.0 >= indices.len() || seeds.1 >= indices.len() {
        return None;
    }
    let project = Projector::new(features, indices);
    match grow(
        &project,
        features,
        indices,
        seeds,
        Trial {
            min_bin_size,
            eligible,
            seed,
        },
    ) {
        Shape::TwoClouds(pieces) => Some(pieces),
        _ => None,
    }
}

fn grow(
    project: &Projector,
    features: &ContigFeatures,
    indices: &[usize],
    seeds: (usize, usize),
    trial: Trial,
) -> Shape {
    let Trial {
        min_bin_size,
        eligible,
        seed,
    } = trial;
    let Some(pieces) = two_means(project, features, indices, seeds) else {
        return Shape::Untestable(Mute::OneSided);
    };
    let first = centroid(features, &pieces[0]);
    let second = centroid(features, &pieces[1]);
    let to_first = project.to(&first);
    let to_second = project.to(&second);
    if !pieces
        .iter()
        .all(|piece| features.bin_size(piece) >= min_bin_size)
    {
        return Shape::Untestable(Mute::PieceTooSmall);
    }
    match bimodal(
        &project.metric,
        &to_first,
        &to_second,
        &first,
        &second,
        eligible,
        seed,
    ) {
        true => Shape::TwoClouds(pieces),
        false => Shape::OneCloud,
    }
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
    let project = Projector::new(features, &indices);
    let first = centroid(features, largest);
    let second = centroid(features, next);
    bimodal(
        &project.metric,
        &project.to(&first),
        &project.to(&second),
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

fn two_means(
    project: &Projector,
    features: &ContigFeatures,
    indices: &[usize],
    seeds: (usize, usize),
) -> Option<[Vec<usize>; 2]> {
    let mut side = project.nearer(seeds.0, seeds.1);
    let mut pieces = members(indices, &side)?;
    for _ in 0..MAX_ROUNDS {
        let first = centroid(features, &pieces[0]);
        let second = centroid(features, &pieces[1]);
        let to_first = project.to(&first);
        let to_second = project.to(&second);
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
    }
    Some(pieces)
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
    if n <= crate::tuning::DIP_SAMPLE {
        return (projection, weights);
    }
    let positions = crate::seeds::sample_positions(n, crate::tuning::DIP_SAMPLE, seed);
    (
        positions.iter().map(|p| projection[*p]).collect(),
        positions.iter().map(|p| weights[*p]).collect(),
    )
}
