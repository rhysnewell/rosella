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
