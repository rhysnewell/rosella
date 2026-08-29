//! The dimensionality estimate exists to replace a rule that keys on the coverage sample
//! count, so what it has to get right is the dimension of the manifold rather than the
//! number of columns the data arrives in.

use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::embedding::knn::build_knn;
use rosella::embedding::metrics::euclidean;

const AMBIENT: usize = 20;

/// Uniform points on a `manifold`-dimensional linear subspace of `AMBIENT` dimensions.
fn embedded_manifold(n_points: usize, manifold: usize, seed: u64) -> Vec<Vec<f64>> {
    let mut rng = StdRng::seed_from_u64(seed);
    let basis: Vec<Vec<f64>> = (0..manifold)
        .map(|_| (0..AMBIENT).map(|_| rng.random_range(-1.0..1.0)).collect())
        .collect();

    (0..n_points)
        .map(|_| {
            let coordinates: Vec<f64> = (0..manifold).map(|_| rng.random_range(0.0..1.0)).collect();
            (0..AMBIENT)
                .map(|column| {
                    coordinates
                        .iter()
                        .zip(basis.iter())
                        .map(|(value, axis)| value * axis[column])
                        .sum()
                })
                .collect()
        })
        .collect()
}

#[test]
fn estimates_the_manifold_dimension_not_the_ambient_one() {
    for manifold in [2usize, 4, 6] {
        let rows = embedded_manifold(3000, manifold, 7);
        let graph = build_knn(rows.len(), 20, 42, |i, j| euclidean(&rows[i], &rows[j]));

        let estimate = graph.intrinsic_dimension().expect("no estimate");
        assert!(
            (estimate - manifold as f64).abs() < 1.0,
            "a {manifold} dimensional manifold in {AMBIENT} dimensions estimated {estimate}"
        );
    }
}

#[test]
fn too_few_points_to_fit_gives_no_estimate() {
    let rows = embedded_manifold(10, 2, 7);
    let graph = build_knn(rows.len(), 5, 42, |i, j| euclidean(&rows[i], &rows[j]));

    assert!(graph.intrinsic_dimension().is_none());
}

#[test]
fn one_neighbour_gives_no_estimate() {
    let rows = embedded_manifold(200, 2, 7);
    let graph = build_knn(rows.len(), 1, 42, |i, j| euclidean(&rows[i], &rows[j]));

    assert!(graph.intrinsic_dimension().is_none());
}
