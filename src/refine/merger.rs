use std::collections::BTreeMap;

use crate::embedding::{features::ContigFeatures, metrics::AggregateMetric};
use crate::refine::bin_stats::{AGGREGATE, bin_stats, centroid};

struct Group {
    indices: Vec<usize>,
    size: usize,
    aggregate: f64,
}

fn union(left: &[usize], right: &[usize]) -> Vec<usize> {
    let mut merged = Vec::with_capacity(left.len() + right.len());
    merged.extend_from_slice(left);
    merged.extend_from_slice(right);
    merged.sort_unstable();
    merged
}

fn bar(left: &Group, right: &Group) -> f64 {
    let total = (left.size + right.size) as f64;
    (left.aggregate * left.size as f64 + right.aggregate * right.size as f64) / total
}

fn root(parent: &mut [usize], mut node: usize) -> usize {
    while parent[node] != node {
        parent[node] = parent[parent[node]];
        node = parent[node];
    }
    node
}

/// The splitter only ever divides, so a genome cut in two by the first clustering had no way
/// back. The bar is each pair's own spread, so a loose pair has to be loose relative to itself.
pub fn merge_bins(
    features: &ContigFeatures,
    bins: BTreeMap<usize, Vec<usize>>,
    max_bin_size: usize,
    seed: u64,
) -> (BTreeMap<usize, Vec<usize>>, usize) {
    let ids = bins.keys().copied().collect::<Vec<_>>();
    let mut groups = Vec::with_capacity(ids.len());
    let mut centroids = Vec::with_capacity(ids.len());

    for id in &ids {
        let indices = &bins[id];
        let Some(stats) = bin_stats(features, indices, seed) else {
            return (bins, 0);
        };
        centroids.push(centroid(features, indices));
        groups.push(Group {
            indices: indices.clone(),
            size: features.bin_size(indices),
            aggregate: stats.mean[AGGREGATE],
        });
    }

    let metric = AggregateMetric::new(features.n_samples() * 2, features.distance_settings());
    let mut candidates = Vec::new();
    for left in 0..groups.len() {
        for right in (left + 1)..groups.len() {
            let distance = metric.distance(
                &centroids[left].row,
                &centroids[right].row,
                centroids[left].floor,
                centroids[right].floor,
            );
            if distance <= bar(&groups[left], &groups[right]) {
                candidates.push((distance, left, right));
            }
        }
    }
    candidates.sort_by(|a, b| a.0.total_cmp(&b.0));

    let mut parent = (0..groups.len()).collect::<Vec<_>>();
    let mut merges = 0;
    for (_, left, right) in candidates {
        let (left, right) = (root(&mut parent, left), root(&mut parent, right));
        if left == right {
            continue;
        }
        let combined = union(&groups[left].indices, &groups[right].indices);
        if features.bin_size(&combined) > max_bin_size {
            continue;
        }
        let Some(stats) = bin_stats(features, &combined, seed) else {
            continue;
        };
        if stats.mean[AGGREGATE] > bar(&groups[left], &groups[right]) {
            continue;
        }

        parent[right] = left;
        groups[left] = Group {
            size: features.bin_size(&combined),
            aggregate: stats.mean[AGGREGATE],
            indices: combined,
        };
        merges += 1;
    }

    if merges == 0 {
        return (bins, 0);
    }

    let mut merged = BTreeMap::new();
    for position in 0..groups.len() {
        if root(&mut parent, position) == position {
            merged.insert(ids[position], std::mem::take(&mut groups[position].indices));
        }
    }
    (merged, merges)
}
