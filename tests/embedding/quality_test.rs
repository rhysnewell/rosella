//! The detector has to separate a layout that kept the neighbourhoods from one that did not,
//! on data with no labels, or it cannot rank an embedding of an assembly nobody has scored.

use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::embedding::knn::build_knn;
use rosella::embedding::metrics::euclidean;
use rosella::embedding::quality::neighbour_preservation;

const AMBIENT: usize = 20;
const POINTS: usize = 1200;
const NEIGHBOURS: usize = 15;

/// A plane laid into 20 dimensions, returned with the coordinates that generated it.
fn plane(seed: u64) -> (Vec<Vec<f64>>, Array2<f64>) {
    let mut rng = StdRng::seed_from_u64(seed);
    let basis: Vec<Vec<f64>> = (0..2)
        .map(|_| (0..AMBIENT).map(|_| rng.random_range(-1.0..1.0)).collect())
        .collect();

    let mut truth = Array2::zeros((POINTS, 2));
    let mut rows = Vec::with_capacity(POINTS);
    for point in 0..POINTS {
        let coordinates = [rng.random_range(0.0..1.0), rng.random_range(0.0..1.0)];
        truth[[point, 0]] = coordinates[0];
        truth[[point, 1]] = coordinates[1];
        rows.push(
            (0..AMBIENT)
                .map(|column| {
                    coordinates
                        .iter()
                        .zip(basis.iter())
                        .map(|(value, axis)| value * axis[column])
                        .sum()
                })
                .collect(),
        );
    }
    (rows, truth)
}

#[test]
fn the_generating_coordinates_score_far_above_a_scrambled_layout() {
    let (rows, truth) = plane(7);
    let source = build_knn(rows.len(), NEIGHBOURS, 42, |i, j| {
        euclidean(&rows[i], &rows[j])
    });

    let mut rng = StdRng::seed_from_u64(11);
    let scrambled =
        Array2::from_shape_fn((POINTS, 2), |_| rng.random_range(0.0..1.0f64));

    let kept = neighbour_preservation(&truth, &source, 42).expect("no score");
    let lost = neighbour_preservation(&scrambled, &source, 42).expect("no score");

    assert!(kept > 0.5, "the true coordinates only kept {kept}");
    assert!(lost < 0.05, "a random layout kept {lost}");
}

#[test]
fn a_layout_that_lost_a_dimension_scores_between_the_two() {
    let (rows, truth) = plane(7);
    let source = build_knn(rows.len(), NEIGHBOURS, 42, |i, j| {
        euclidean(&rows[i], &rows[j])
    });

    let mut flattened = truth.clone();
    flattened.column_mut(1).fill(0.0);

    let kept = neighbour_preservation(&truth, &source, 42).expect("no score");
    let squashed = neighbour_preservation(&flattened, &source, 42).expect("no score");

    assert!(
        squashed < kept,
        "collapsing an axis scored {squashed} against {kept}"
    );
}
