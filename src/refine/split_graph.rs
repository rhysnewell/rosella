use log::debug;

use crate::embedding::features::ContigFeatures;
use crate::embedding::umap::EmbedOverrides;
use crate::embedding::{Graph, KNN_SPLIT, induced};
use crate::seeds::Seeds;

/// A slice of the whole-assembly graph only serves a bin when it is at least as connected as the
/// graph the bin would have built for itself, so the bar is what `knn_size` would have asked for.
pub fn bin_graph(
    features: &ContigFeatures,
    assembly: Option<&Graph>,
    indices: &[usize],
    n_neighbours: usize,
    seeds: Seeds,
    overrides: &EmbedOverrides,
) -> Graph {
    let own = || features.graph_of(indices, n_neighbours, seeds, overrides, KNN_SPLIT);
    let Some(assembly) = assembly else {
        return own();
    };
    let sliced = induced(assembly, indices);
    if sliced.nnz() < features.knn_size(indices.len(), n_neighbours) * indices.len() {
        debug!("Induced graph too sparse for {} contigs", indices.len());
        return own();
    }
    sliced
}
