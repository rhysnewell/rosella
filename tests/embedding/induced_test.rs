use ndarray::Array2;
use rosella::embedding::knn::KnnGraph;
use rosella::embedding::{Graph, induced};
use sprs::TriMatI;

fn graph(n: usize, edges: &[(usize, usize, f32)]) -> Graph {
    let mut triplets = TriMatI::<f32, u32>::new((n, n));
    for (a, b, weight) in edges {
        triplets.add_triplet(*a, *b, *weight);
        triplets.add_triplet(*b, *a, *weight);
    }
    triplets.to_csr()
}

#[test]
fn induction_renumbers_and_drops_edges_leaving_the_set() {
    let whole = graph(6, &[(1, 3, 0.5), (3, 5, 0.25), (1, 2, 0.75), (0, 4, 1.0)]);

    let part = induced(&whole, &[1, 3, 5]);

    assert_eq!(part.rows(), 3);
    let dense = part.to_dense();
    assert_eq!(dense[[0, 1]], 0.5);
    assert_eq!(dense[[1, 2]], 0.25);
    assert_eq!(dense.iter().filter(|weight| **weight > 0.0).count(), 4);
}

#[test]
fn neighbour_induction_takes_the_widest_width_every_row_can_fill() {
    let indices =
        Array2::from_shape_vec((4, 3), vec![1, 2, 3, 0, 2, 3, 3, 0, 1, 2, 1, 0]).expect("shape");
    let dists = Array2::from_shape_vec(
        (4, 3),
        vec![0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3],
    )
    .expect("shape");
    let knn = KnnGraph { indices, dists };

    let part = knn.induced(&[1, 2, 3]).expect("two survivors a row");

    assert_eq!(part.indices.ncols(), 2);
    assert_eq!(part.indices.row(0).to_vec(), vec![1, 2]);
    assert_eq!(part.dists.row(0).to_vec(), vec![0.2, 0.3]);
    assert_eq!(part.indices.row(2).to_vec(), vec![1, 0]);
}

#[test]
fn neighbour_induction_refuses_a_row_the_set_leaves_alone() {
    let indices = Array2::from_shape_vec((4, 2), vec![1, 2, 0, 3, 3, 0, 1, 2]).expect("shape");
    let dists = Array2::from_shape_vec((4, 2), vec![0.1, 0.2, 0.1, 0.4, 0.3, 0.4, 0.2, 0.3])
        .expect("shape");
    let knn = KnnGraph { indices, dists };

    assert!(knn.induced(&[0, 1, 2, 3]).is_some());
    assert!(knn.induced(&[0, 1]).is_none());
}
