use ndarray::{ArrayBase, Data, Ix2};
use rayon::prelude::*;

const NO_CLUSTERS: f64 = -1.0;
const EPSILON: f64 = 1e-12;

/// Density Based Cluster Validity (Moulavi et al. 2014), the score flight ranks its
/// HDBSCAN parameter sweep by. Ranges from -1 to 1, higher is better.
pub fn dbcv<S: Data<Elem = f64> + Sync>(points: &ArrayBase<S, Ix2>, labels: &[i32]) -> f64 {
    let dimensionality = points.ncols() as f64;
    let clusters = group_by_label(labels);
    if clusters.len() < 2 {
        return NO_CLUSTERS;
    }

    let geometries = clusters
        .par_iter()
        .map(|indices| ClusterGeometry::new(points, indices, dimensionality))
        .collect::<Vec<_>>();

    let scores = (0..geometries.len())
        .into_par_iter()
        .map(|i| {
            let separation = (0..geometries.len())
                .filter(|j| *j != i)
                .map(|j| density_separation(points, &geometries[i], &geometries[j]))
                .fold(f64::INFINITY, f64::min);

            let sparseness = geometries[i].sparseness;
            let denominator = separation.max(sparseness);
            let validity = if denominator <= EPSILON {
                0.0
            } else {
                (separation - sparseness) / denominator
            };

            geometries[i].indices.len() as f64 / labels.len() as f64 * validity
        })
        .sum::<f64>();

    if scores.is_nan() { NO_CLUSTERS } else { scores }
}

fn group_by_label(labels: &[i32]) -> Vec<Vec<usize>> {
    let highest = labels.iter().copied().max().unwrap_or(-1);
    if highest < 0 {
        return Vec::new();
    }

    let mut clusters = vec![Vec::new(); highest as usize + 1];
    for (index, label) in labels.iter().enumerate() {
        if *label >= 0 {
            clusters[*label as usize].push(index);
        }
    }
    clusters.retain(|cluster| cluster.len() > 1);
    clusters
}

struct ClusterGeometry {
    indices: Vec<usize>,
    core_distances: Vec<f64>,
    /// Nodes with more than one edge in the cluster's minimum spanning tree. Leaves are
    /// excluded so a single stray point cannot dominate the score.
    internal: Vec<bool>,
    sparseness: f64,
}

impl ClusterGeometry {
    fn new<S: Data<Elem = f64> + Sync>(
        points: &ArrayBase<S, Ix2>,
        indices: &[usize],
        dimensionality: f64,
    ) -> Self {
        let core_distances = all_points_core_distances(points, indices, dimensionality);
        let edges = minimum_spanning_tree(points, indices, &core_distances);

        let mut degree = vec![0usize; indices.len()];
        for (from, to, _) in edges.iter() {
            degree[*from] += 1;
            degree[*to] += 1;
        }
        let internal = degree.iter().map(|d| *d > 1).collect::<Vec<_>>();

        let internal_max = edges
            .iter()
            .filter(|(from, to, _)| internal[*from] && internal[*to])
            .map(|(_, _, weight)| *weight)
            .fold(f64::NEG_INFINITY, f64::max);
        let overall_max = edges
            .iter()
            .map(|(_, _, weight)| *weight)
            .fold(f64::NEG_INFINITY, f64::max);

        let sparseness = if internal_max.is_finite() {
            internal_max
        } else if overall_max.is_finite() {
            overall_max
        } else {
            0.0
        };

        Self {
            indices: indices.to_vec(),
            core_distances,
            internal,
            sparseness,
        }
    }
}

/// The all points core distance: the inverse of the mean inverse distance to every other
/// member of the cluster, raised to the dimensionality.
fn all_points_core_distances<S: Data<Elem = f64> + Sync>(
    points: &ArrayBase<S, Ix2>,
    indices: &[usize],
    dimensionality: f64,
) -> Vec<f64> {
    indices
        .par_iter()
        .map(|i| {
            let total = indices
                .iter()
                .filter(|j| *j != i)
                .map(|j| {
                    let distance = euclidean(points, *i, *j).max(EPSILON);
                    distance.powf(-dimensionality)
                })
                .sum::<f64>();

            let mean = total / (indices.len() - 1) as f64;
            if mean <= 0.0 {
                0.0
            } else {
                mean.powf(-1.0 / dimensionality)
            }
        })
        .collect()
}

/// Prim's algorithm over the mutual reachability graph, computing weights as it goes so
/// no square distance matrix is ever held.
///
/// Mutual reachability flattens every short edge to the larger core distance, so ties are
/// the rule rather than the exception and the tree is not unique. Ties break on the raw
/// distance, which both fixes the tree and connects each point to its nearest equally
/// reachable neighbour.
fn minimum_spanning_tree<S: Data<Elem = f64> + Sync>(
    points: &ArrayBase<S, Ix2>,
    indices: &[usize],
    core_distances: &[f64],
) -> Vec<(usize, usize, f64)> {
    let n = indices.len();
    let mut in_tree = vec![false; n];
    let mut best_weight = vec![f64::INFINITY; n];
    let mut best_distance = vec![f64::INFINITY; n];
    let mut best_source = vec![0usize; n];
    let mut edges = Vec::with_capacity(n.saturating_sub(1));

    best_weight[0] = 0.0;
    best_distance[0] = 0.0;
    for step in 0..n {
        let next = (0..n)
            .filter(|position| !in_tree[*position])
            .min_by(|a, b| {
                best_weight[*a]
                    .total_cmp(&best_weight[*b])
                    .then(best_distance[*a].total_cmp(&best_distance[*b]))
            });
        let Some(next) = next else { break };

        in_tree[next] = true;
        if step > 0 {
            edges.push((best_source[next], next, best_weight[next]));
        }

        for other in 0..n {
            if in_tree[other] {
                continue;
            }
            let distance = euclidean(points, indices[next], indices[other]);
            let weight = distance.max(core_distances[next]).max(core_distances[other]);
            let improves = weight < best_weight[other]
                || (weight == best_weight[other] && distance < best_distance[other]);
            if improves {
                best_weight[other] = weight;
                best_distance[other] = distance;
                best_source[other] = next;
            }
        }
    }

    edges
}

/// Closest mutual reachability between two clusters, over their internal nodes only.
fn density_separation<S: Data<Elem = f64> + Sync>(
    points: &ArrayBase<S, Ix2>,
    left: &ClusterGeometry,
    right: &ClusterGeometry,
) -> f64 {
    let candidates = |geometry: &ClusterGeometry| {
        let internal = (0..geometry.indices.len())
            .filter(|position| geometry.internal[*position])
            .collect::<Vec<_>>();
        if internal.is_empty() {
            (0..geometry.indices.len()).collect()
        } else {
            internal
        }
    };

    let left_positions = candidates(left);
    let right_positions = candidates(right);

    left_positions
        .par_iter()
        .map(|a| {
            right_positions
                .iter()
                .map(|b| {
                    euclidean(points, left.indices[*a], right.indices[*b])
                        .max(left.core_distances[*a])
                        .max(right.core_distances[*b])
                })
                .fold(f64::INFINITY, f64::min)
        })
        .reduce(|| f64::INFINITY, f64::min)
}

fn euclidean<S: Data<Elem = f64> + Sync>(points: &ArrayBase<S, Ix2>, a: usize, b: usize) -> f64 {
    points
        .row(a)
        .iter()
        .zip(points.row(b).iter())
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt()
}
