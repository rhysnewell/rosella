use std::collections::BTreeMap;

use crate::embedding::{features::ContigFeatures, metrics::AggregateMetric};
use crate::refine::bin_stats::{AGGREGATE, BinStats, bin_stats, centroid};
use crate::refine::solo;

pub const MERGE_BAR_NAMES: [&str; 2] = ["pair", "widest"];

/// `Widest` reads the loosest contig a bin already holds, so the scale comes from the contigs
/// in front of it rather than one number the whole run is held to.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum MergeBar {
    #[default]
    Pair,
    Widest,
}

impl MergeBar {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "pair" => Some(Self::Pair),
            "widest" => Some(Self::Widest),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct MergeSettings {
    pub bar: MergeBar,
    pub mutual: bool,
    pub short_side: bool,
    pub genome_floor: Option<usize>,
    pub max_bin_size: usize,
    pub seed: u64,
}

struct Group {
    indices: Vec<usize>,
    size: usize,
    mean: Option<f64>,
    widest: Option<f64>,
}

impl Group {
    fn new(features: &ContigFeatures, indices: Vec<usize>, stats: Option<BinStats>) -> Self {
        Self {
            size: features.bin_size(&indices),
            mean: stats.as_ref().map(|stats| stats.mean[AGGREGATE]),
            widest: stats.as_ref().and_then(|stats| {
                stats
                    .per_contig
                    .iter()
                    .map(|row| row[AGGREGATE])
                    .max_by(f64::total_cmp)
            }),
            indices,
        }
    }

    fn spread(&self, bar: MergeBar) -> Option<f64> {
        match bar {
            MergeBar::Pair => self.mean,
            MergeBar::Widest => self.widest,
        }
    }
}

fn union(left: &[usize], right: &[usize]) -> Vec<usize> {
    let mut merged = Vec::with_capacity(left.len() + right.len());
    merged.extend_from_slice(left);
    merged.extend_from_slice(right);
    merged.sort_unstable();
    merged
}

/// A bin of one contig has no spread of its own, so it is offered the one its partner already
/// tolerates. Two of them have nothing to be judged against.
fn bar(left: &Group, right: &Group, mode: MergeBar) -> Option<f64> {
    match (left.spread(mode), right.spread(mode)) {
        (Some(one), Some(other)) => Some(match mode {
            MergeBar::Pair => {
                let total = (left.size + right.size) as f64;
                (one * left.size as f64 + other * right.size as f64) / total
            }
            MergeBar::Widest => one.max(other),
        }),
        (Some(only), None) | (None, Some(only)) => Some(only),
        (None, None) => None,
    }
}

/// A bin already holding a whole genome has nothing to gain from a partner, so only a short
/// side is allowed to recruit.
fn short_side(left: &Group, right: &Group, genome: usize) -> bool {
    left.size < genome || right.size < genome
}

fn reciprocated(best: &[Option<(f64, usize)>]) -> Vec<(f64, usize, usize)> {
    let mut pairs = Vec::new();
    for (left, nearest) in best.iter().enumerate() {
        let Some((distance, right)) = *nearest else {
            continue;
        };
        if left < right && best[right].is_some_and(|(_, back)| back == left) {
            pairs.push((distance, left, right));
        }
    }
    pairs
}

fn offer(best: &mut [Option<(f64, usize)>], one: usize, other: usize, distance: f64) {
    if best[one].is_none_or(|(closest, _)| distance < closest) {
        best[one] = Some((distance, other));
    }
}

/// The splitter only ever divides, so a genome cut in two by the first clustering had no way
/// back. A genome floor admits bins of one contig, which have no spread, and refuses any union
/// solo would take apart again.
pub fn merge_bins(
    features: &ContigFeatures,
    bins: BTreeMap<usize, Vec<usize>>,
    settings: MergeSettings,
) -> (BTreeMap<usize, Vec<usize>>, usize) {
    let mut ids = Vec::with_capacity(bins.len());
    let mut groups = Vec::with_capacity(bins.len());
    let mut centroids = Vec::with_capacity(bins.len());
    let mut unscored = BTreeMap::new();

    for (id, indices) in bins {
        let stats = bin_stats(features, &indices, settings.seed);
        if stats.is_none() && (settings.genome_floor.is_none() || indices.is_empty()) {
            unscored.insert(id, indices);
            continue;
        }
        ids.push(id);
        centroids.push(centroid(features, &indices));
        groups.push(Group::new(features, indices, stats));
    }

    let genome = settings.genome_floor.map(|floor| floor * 2);
    let metric = AggregateMetric::new(features.n_samples() * 2, features.distance_settings())
        .with_bands(features.bands());
    let mut best = vec![None; groups.len()];
    let mut candidates = Vec::new();
    for left in 0..groups.len() {
        for right in (left + 1)..groups.len() {
            let Some(bar) = bar(&groups[left], &groups[right], settings.bar) else {
                continue;
            };
            if settings.short_side
                && genome.is_some_and(|genome| !short_side(&groups[left], &groups[right], genome))
            {
                continue;
            }
            let distance = metric.distance(
                &centroids[left].row,
                &centroids[right].row,
                centroids[left].floor,
                centroids[right].floor,
            );
            if settings.mutual {
                offer(&mut best, left, right, distance);
                offer(&mut best, right, left, distance);
            } else if distance <= bar {
                candidates.push((distance, left, right));
            }
        }
    }
    if settings.mutual {
        candidates = reciprocated(&best);
    }
    candidates.sort_by(|a, b| a.0.total_cmp(&b.0));

    let mut parent = (0..groups.len()).collect::<Vec<_>>();
    let mut merges = 0;
    for (_, left, right) in candidates {
        let (left, right) = (root(&mut parent, left), root(&mut parent, right));
        if left == right {
            continue;
        }
        let Some(bar) = bar(&groups[left], &groups[right], settings.bar) else {
            continue;
        };
        let combined = union(&groups[left].indices, &groups[right].indices);
        if features.bin_size(&combined) > settings.max_bin_size {
            continue;
        }
        if settings
            .genome_floor
            .is_some_and(|floor| solo::candidate(features, &combined, floor, false).is_some())
        {
            continue;
        }
        let Some(stats) = bin_stats(features, &combined, settings.seed) else {
            continue;
        };
        if stats.mean[AGGREGATE] > bar {
            continue;
        }

        parent[right] = left;
        groups[left] = Group::new(features, combined, Some(stats));
        merges += 1;
    }

    let mut merged = unscored;
    for position in 0..groups.len() {
        if root(&mut parent, position) == position {
            merged.insert(ids[position], std::mem::take(&mut groups[position].indices));
        }
    }
    (merged, merges)
}

fn root(parent: &mut [usize], mut node: usize) -> usize {
    while parent[node] != node {
        parent[node] = parent[parent[node]];
        node = parent[node];
    }
    node
}
