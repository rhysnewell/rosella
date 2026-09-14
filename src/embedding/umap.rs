use ndarray::Array2;
use umap_rs::{GraphParams, ManifoldParams, Umap, UmapConfig};

use crate::embedding::{Graph, knn::KnnGraph};

/// Set from the CLI so an ablation can hold the neighbour search still.
#[derive(Debug, Clone, Copy, Default)]
pub struct EmbedOverrides {
    pub knn_candidates: Option<usize>,
    pub graph_weights: crate::embedding::manifold::GraphWeights,
}

pub fn manifold_graph(n_points: usize, knn: &KnnGraph, n_neighbours: usize) -> Graph {
    let config = UmapConfig {
        manifold: ManifoldParams::default(),
        graph: GraphParams {
            n_neighbors: n_neighbours.min(knn.indices.ncols()),
            set_op_mix_ratio: 1.0,
            ..Default::default()
        },
        ..Default::default()
    };

    let _timer = crate::timing::scope("manifold");
    let empty = Array2::<f32>::zeros((n_points, 0));
    let manifold =
        Umap::new(config).learn_manifold(empty.view(), knn.indices.view(), knn.dists.view());
    manifold.graph().clone()
}
