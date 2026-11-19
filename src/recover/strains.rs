use std::collections::BTreeMap;

use ndarray::ArrayView2;

use crate::embedding::features::ContigFeatures;
use crate::embedding::metrics::AggregateMetric;
use crate::markers::replicon::DEPTH_FLOOR;
use crate::quality::Scorer;
use crate::recover::coherence::Coherence;
use crate::recover::fold::HOST;
use crate::refine::bin_stats::centroid;

const ROUNDS: usize = 32;

pub struct Depths<'a> {
    pub depths: ArrayView2<'a, f64>,
    pub lengths: &'a [usize],
    pub read_lengths: Option<&'a [f64]>,
}

impl Depths<'_> {
    fn logs(&self, contig: usize) -> Vec<f64> {
        self.depths
            .row(contig)
            .iter()
            .map(|depth| (depth + DEPTH_FLOOR).ln())
            .collect()
    }

    fn centre(&self, contigs: &[usize]) -> Vec<f64> {
        let total = contigs
            .iter()
            .map(|contig| self.lengths[*contig])
            .sum::<usize>() as f64;
        let mut mean = vec![0.0; self.depths.ncols()];
        for contig in contigs {
            let weight = self.lengths[*contig] as f64 / total;
            for (sum, log) in mean.iter_mut().zip(self.logs(*contig)) {
                *sum += weight * log;
            }
        }
        mean
    }

    fn fit(&self, contig: usize, centre: &[f64], coherence: &Coherence) -> f64 {
        let length = self.lengths[contig] as f64;
        let spread = coherence.spread(self.lengths[contig]);
        self.logs(contig)
            .iter()
            .zip(centre)
            .zip(spread)
            .enumerate()
            .filter(|(_, (_, spread))| **spread > 0.0)
            .map(|(sample, ((log, centre), spread))| {
                let depth = centre.exp() - DEPTH_FLOOR;
                let count = self.read_lengths.map_or(0.0, |reads| {
                    depth.max(0.0) * reads[sample] / length / (depth.max(0.0) + DEPTH_FLOOR).powi(2)
                });
                let variance = spread + count;
                -0.5 * ((log - centre).powi(2) / variance + variance.ln())
            })
            .sum()
    }
}

fn leading(differences: &[Vec<f64>]) -> Vec<f64> {
    let width = differences[0].len();
    let mut axis = differences
        .iter()
        .max_by(|a, b| norm(a).total_cmp(&norm(b)))
        .cloned()
        .unwrap_or_else(|| vec![1.0; width]);
    for _ in 0..ROUNDS {
        let mut next = vec![0.0; width];
        for difference in differences {
            let along = dot(difference, &axis);
            for (sum, value) in next.iter_mut().zip(difference) {
                *sum += along * value;
            }
        }
        let size = norm(&next);
        if size == 0.0 {
            break;
        }
        axis = next.iter().map(|value| value / size).collect();
    }
    axis
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

// Two strains carry one copy each of a single-copy marker, so a doubled marker's two carriers
// say which depth each strain sits at, and the markers only judge the cut depth made.
fn halves(
    members: &[usize],
    carriers: &[Vec<u32>],
    depths: &Depths,
    coherence: &Coherence,
    features: &ContigFeatures,
) -> Option<(Vec<usize>, Vec<usize>)> {
    let pairs = carriers
        .iter()
        .filter(|family| family.len() == 2)
        .map(|family| (members[family[0] as usize], members[family[1] as usize]))
        .collect::<Vec<_>>();
    if pairs.is_empty() {
        return None;
    }
    let (high, low): (Vec<_>, Vec<_>) = oriented(&pairs, depths, coherence).into_iter().unzip();
    let (upper, lower) = assigned(members, &high, &low, features);
    (!upper.is_empty() && !lower.is_empty()).then_some((upper, lower))
}

// A low-depth sample's read count noise can outweigh a strain's depth offset, so the axis only
// starts the orientation and each pair then settles on the noise-weighted fit.
fn oriented(
    pairs: &[(usize, usize)],
    depths: &Depths,
    coherence: &Coherence,
) -> Vec<(usize, usize)> {
    let differences = pairs
        .iter()
        .map(|(a, b)| {
            depths
                .logs(*a)
                .iter()
                .zip(depths.logs(*b))
                .map(|(x, y)| x - y)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let axis = leading(&differences);
    let mut pairs = pairs
        .iter()
        .zip(&differences)
        .map(|((a, b), difference)| match dot(difference, &axis) >= 0.0 {
            true => (*a, *b),
            false => (*b, *a),
        })
        .collect::<Vec<_>>();
    for _ in 0..ROUNDS {
        let (high, low): (Vec<_>, Vec<_>) = pairs.iter().copied().unzip();
        let (high, low) = (depths.centre(&high), depths.centre(&low));
        let fit = |contig, centre: &[f64]| depths.fit(contig, centre, coherence);
        let next = pairs
            .iter()
            .map(
                |(a, b)| match fit(*a, &high) + fit(*b, &low) >= fit(*b, &high) + fit(*a, &low) {
                    true => (*a, *b),
                    false => (*b, *a),
                },
            )
            .collect::<Vec<_>>();
        if next == pairs {
            break;
        }
        pairs = next;
    }
    pairs
}

fn assigned(
    members: &[usize],
    high: &[usize],
    low: &[usize],
    features: &ContigFeatures,
) -> (Vec<usize>, Vec<usize>) {
    let metric = AggregateMetric::new(
        features.n_samples() * 2,
        features.distance_settings().clone(),
    );
    let points = features.points(members);
    let (mut upper, mut lower) = (high.to_vec(), low.to_vec());
    for _ in 0..ROUNDS {
        let (first, second) = (centroid(features, &upper), centroid(features, &lower));
        let (mut next_upper, mut next_lower) = (Vec::new(), Vec::new());
        for (contig, point) in members.iter().zip(&points) {
            match metric.distance(point, &first) <= metric.distance(point, &second) {
                true => next_upper.push(*contig),
                false => next_lower.push(*contig),
            }
        }
        let settled = next_upper == upper;
        (upper, lower) = (next_upper, next_lower);
        if settled || upper.is_empty() || lower.is_empty() {
            break;
        }
    }
    (upper, lower)
}

// One genome cut in two cannot hold more markers than its whole. Two strains hold a medium
// quality genome's worth more.
pub fn split_strains(
    bins: &mut BTreeMap<usize, Vec<usize>>,
    quality: &dyn Scorer,
    carriers: impl Fn(&[usize]) -> Vec<Vec<u32>>,
    depths: &Depths,
    coherence: &Coherence,
    features: &ContigFeatures,
) -> usize {
    let mut next = bins.keys().max().map_or(0, |last| last + 1);
    let mut queue = bins.keys().copied().collect::<Vec<_>>();
    let mut splits = 0;
    while let Some(bin) = queue.pop() {
        let members = &bins[&bin];
        let whole = quality.score(members);
        if whole.contamination <= HOST.contamination {
            continue;
        }
        let Some((upper, lower)) = halves(members, &carriers(members), depths, coherence, features)
        else {
            continue;
        };
        let (first, second) = (quality.score(&upper), quality.score(&lower));
        let gained = first.completeness + second.completeness - whole.completeness;
        if first.completeness < HOST.completeness
            || second.completeness < HOST.completeness
            || gained < HOST.completeness
        {
            continue;
        }
        bins.insert(bin, upper);
        bins.insert(next, lower);
        queue.extend([bin, next]);
        next += 1;
        splits += 1;
    }
    splits
}
