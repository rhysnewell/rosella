pub mod features;
pub mod knn;
pub mod manifold;
pub mod metrics;
pub mod umap;

pub const KNN_ASSEMBLY: &str = "knn_assembly";
pub const KNN_POOL: &str = "knn_pool";
pub const KNN_SPLIT: &str = "knn_split";

pub type Graph = sprs::CsMatI<f32, u32, usize>;

pub(crate) fn row_of(graph: &Graph, row: usize) -> (&[u32], &[f32]) {
    let start = graph.indptr().index(row);
    let end = graph.indptr().index(row + 1);
    (&graph.indices()[start..end], &graph.data()[start..end])
}

/// `Level::from_graph` and `label_propagation` size their vectors by `graph.rows()` and index
/// them by node id, so a slice has to come back renumbered rather than masked.
pub fn induced(graph: &Graph, indices: &[usize]) -> Graph {
    let mut triplets = sprs::TriMatI::<f32, u32>::new((indices.len(), indices.len()));
    for (row, node) in indices.iter().enumerate() {
        let (columns, weights) = row_of(graph, *node);
        for (column, weight) in columns.iter().zip(weights) {
            if let Ok(position) = indices.binary_search(&(*column as usize)) {
                triplets.add_triplet(row, position, *weight);
            }
        }
    }
    triplets.to_csr()
}
