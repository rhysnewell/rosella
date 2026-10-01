use crate::embedding::{
    features::ContigFeatures,
    metrics::{AggregateMetric, Point},
};
use crate::refine::bar::MIN_SPLIT_CONTIGS;
use crate::refine::bin_stats::centroid;
use crate::refine::dip;

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
    points: Vec<Point>,
}

impl Projector {
    fn new(features: &ContigFeatures, indices: &[usize]) -> Self {
        Self {
            metric: AggregateMetric::new(features.n_samples() * 2, features.distance_settings()),
            points: features.points(indices),
        }
    }

    fn to(&self, centre: &Point) -> Vec<f64> {
        self.points
            .iter()
            .map(|point| self.metric.distance(point, centre))
            .collect()
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
    first: &Point,
    second: &Point,
    eligible: usize,
    seed: u64,
) -> bool {
    // The difference of two distances saturates at the centroid gap past either centroid,
    // so it piles any cloud up at both ends. The coordinate along the axis does not.
    let gap = metric.distance(first, second);
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
    dip::exceeds_null(&tested(projection, seed), draws_for(eligible), seed)
}

fn tested(projection: Vec<f64>, seed: u64) -> Vec<f64> {
    let n = projection.len();
    if n <= crate::tuning::DIP_SAMPLE {
        return projection;
    }
    crate::seeds::sample_positions(n, crate::tuning::DIP_SAMPLE, seed)
        .iter()
        .map(|position| projection[*position])
        .collect()
}
