use ndarray::array;
use rosella::embedding::fuzzy::{Kernel, reshaped_graph};
use rosella::embedding::knn::KnnGraph;

#[test]
fn a_global_kernel_weighs_by_distance_and_drops_edges_past_the_cutoff() {
    let knn = KnnGraph {
        indices: array![[1, 2], [0, 2], [0, 1]],
        dists: array![[0.1, 0.4], [0.1, 0.3], [0.4, 0.3]],
    };

    let dense = reshaped_graph(3, &knn, Kernel::Global(0.1), 0.3).to_dense();

    let near = (-1.0f32).exp();
    let far = (-3.0f32).exp();
    assert!((dense[[0, 1]] - (2.0 * near - near * near)).abs() < 1e-6);
    assert!((dense[[1, 2]] - (2.0 * far - far * far)).abs() < 1e-6);
    assert_eq!(dense[[0, 2]], 0.0);
    assert_eq!(dense[[2, 0]], 0.0);
}
