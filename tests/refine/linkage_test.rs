use ndarray::Array2;
use rosella::embedding::knn::KnnGraph;
use rosella::refine::linkage::candidates;

const PIECE: usize = 100_000;

fn chain(gaps: &[f32]) -> KnnGraph {
    let points = gaps.len() + 1;
    let mut indices = Array2::<u32>::zeros((points, 2));
    let mut dists = Array2::<f32>::from_elem((points, 2), f32::MAX);
    for node in 0..points {
        indices[[node, 0]] = node as u32;
        dists[[node, 0]] = 0.0;
        if node + 1 < points {
            indices[[node, 1]] = node as u32 + 1;
            dists[[node, 1]] = gaps[node];
        } else {
            indices[[node, 1]] = node as u32;
        }
    }
    KnnGraph { indices, dists }
}

fn order(points: usize) -> Vec<usize> {
    (0..points).collect()
}

#[test]
fn the_tight_group_is_offered_before_the_rest_joins_it() {
    let knn = chain(&[0.01, 0.01, 0.01, 0.9, 0.02, 0.02]);
    let found = candidates(&knn, &order(7), |_| PIECE, 2 * PIECE, 10 * PIECE);

    assert!(found.contains(&vec![0, 1, 2, 3]));
    assert!(found.contains(&vec![4, 5, 6]));
    assert!(found.contains(&vec![0, 1, 2, 3, 4, 5, 6]));
}

/// A component under the run's genome floor is a shard, and one over the bin ceiling is the
/// whole pool, so neither is a proposal worth a prediction.
#[test]
fn the_floor_and_the_ceiling_both_bind() {
    let knn = chain(&[0.01, 0.02, 0.03]);
    let found = candidates(&knn, &order(4), |_| PIECE, 3 * PIECE, 3 * PIECE);

    assert_eq!(found, vec![vec![0, 1, 2]]);
}
