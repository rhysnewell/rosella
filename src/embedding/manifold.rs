use sprs::TriMatI;

use crate::embedding::{Graph, knn::KnnGraph};

pub const GRAPH_WEIGHT_NAMES: [&str; 3] = ["fuzzy", "snn", "local-scale"];

/// How a kNN graph becomes the weighted graph a partition reads. `Fuzzy` is UMAP's smooth
/// kNN sigma search and fuzzy union; the other two read the same neighbour lists directly.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum GraphWeights {
    #[default]
    Fuzzy,
    Snn,
    LocalScale,
}

impl GraphWeights {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "fuzzy" => Some(Self::Fuzzy),
            "snn" => Some(Self::Snn),
            "local-scale" => Some(Self::LocalScale),
            _ => None,
        }
    }
}

/// Jarvis and Patrick (1973). Two contigs are close when they keep the same company, which a
/// scalar distance cannot say.
pub fn shared_neighbours(knn: &KnnGraph) -> Graph {
    let _timer = crate::timing::scope("manifold");
    let sorted = sorted_neighbours(knn);
    build(knn, &sorted, |_, a, b| jaccard(&sorted[a], &sorted[b]))
}

/// Zelnik-Manor and Perona (2004). The distance to the k-th neighbour stands in for the local
/// density, so a tight cluster and a diffuse one are read on their own scales.
pub fn local_scaled(knn: &KnnGraph) -> Graph {
    let _timer = crate::timing::scope("manifold");
    let sorted = sorted_neighbours(knn);
    let scales = scales(knn);
    build(knn, &sorted, |distance, a, b| {
        let scale = scales[a] * scales[b];
        if scale <= 0.0 {
            0.0
        } else {
            (-(distance * distance) / scale).exp()
        }
    })
}

/// The descent returns each row ordered by distance, and both weights ask set questions of
/// the neighbour lists, so they get one index-ordered copy to binary search.
fn sorted_neighbours(knn: &KnnGraph) -> Vec<Vec<u32>> {
    (0..knn.n_points())
        .map(|row| {
            let mut neighbours = knn
                .indices
                .row(row)
                .iter()
                .copied()
                .filter(|column| *column as usize != row)
                .collect::<Vec<_>>();
            neighbours.sort_unstable();
            neighbours.dedup();
            neighbours
        })
        .collect()
}

/// Only mutual edges survive, which is what makes both of these cheaper than the fuzzy union:
/// a neighbour one end does not return is a neighbour neither end believes in.
fn build<W>(knn: &KnnGraph, sorted: &[Vec<u32>], weight: W) -> Graph
where
    W: Fn(f32, usize, usize) -> f32,
{
    let n_points = knn.n_points();
    let mut triplets = TriMatI::<f32, u32>::new((n_points, n_points));

    for row in 0..n_points {
        for (column, distance) in knn
            .indices
            .row(row)
            .iter()
            .zip(knn.dists.row(row).iter())
            .map(|(column, distance)| (*column as usize, *distance))
        {
            if column <= row || sorted[column].binary_search(&(row as u32)).is_err() {
                continue;
            }
            let value = weight(distance, row, column);
            if value > 0.0 {
                triplets.add_triplet(row, column, value);
                triplets.add_triplet(column, row, value);
            }
        }
    }

    triplets.to_csr()
}

fn jaccard(a: &[u32], b: &[u32]) -> f32 {
    let mut shared = 0usize;
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                shared += 1;
                i += 1;
                j += 1;
            }
        }
    }
    let union = a.len() + b.len() - shared;
    if union == 0 {
        0.0
    } else {
        shared as f32 / union as f32
    }
}

fn scales(knn: &KnnGraph) -> Vec<f32> {
    let last = knn.indices.ncols().saturating_sub(1);
    (0..knn.n_points())
        .map(|row| {
            let scale = knn.dists[[row, last]];
            if scale.is_finite() && scale > 0.0 {
                scale
            } else {
                0.0
            }
        })
        .collect()
}
