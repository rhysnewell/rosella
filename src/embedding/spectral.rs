use ndarray::{Array1, Array2};
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::embedding::Graph;

const POWER_ITERATIONS: usize = 200;
const TARGET_SPAN: f32 = 10.0;
const JITTER: f32 = 1e-4;

/// Eigenvectors of the normalised graph Laplacian, which is what UMAP initialises from.
/// A uniform random start leaves flight's low `b` clumping the embedding into hundreds of
/// disconnected knots that the optimiser never merges.
pub fn spectral_init(graph: &Graph, n_components: usize, seed: u64) -> Array2<f32> {
    let n = graph.rows();
    if n == 0 {
        return Array2::zeros((0, n_components));
    }

    let inverse_sqrt_degree = degrees(graph).mapv(|d| if d > 0.0 { 1.0 / d.sqrt() } else { 0.0 });

    // The all-ones vector scaled by the degrees is the Laplacian's trivial eigenvector.
    // It carries no layout information, so it gets projected out every iteration.
    let mut trivial = degrees(graph).mapv(|d| d.sqrt());
    normalise(&mut trivial);

    let mut basis = random_basis(n, n_components, seed);
    for _ in 0..POWER_ITERATIONS {
        basis = shifted_multiply(graph, &inverse_sqrt_degree, &basis);
        deflate(&mut basis, &trivial);
        orthonormalise(&mut basis);
    }

    rescale(&mut basis, seed);
    basis
}

fn degrees(graph: &Graph) -> Array1<f32> {
    let mut degrees = Array1::zeros(graph.rows());
    for row in 0..graph.rows() {
        let start = graph.indptr().index(row);
        let end = graph.indptr().index(row + 1);
        degrees[row] = graph.data()[start..end].iter().sum();
    }
    degrees
}

/// `(I + D^-1/2 A D^-1/2) X`. The shift keeps every eigenvalue non-negative so power
/// iteration converges on the largest rather than the most negative.
fn shifted_multiply(
    graph: &Graph,
    inverse_sqrt_degree: &Array1<f32>,
    basis: &Array2<f32>,
) -> Array2<f32> {
    let n_components = basis.ncols();
    let rows = (0..graph.rows())
        .into_par_iter()
        .map(|row| {
            let start = graph.indptr().index(row);
            let end = graph.indptr().index(row + 1);
            let scale = inverse_sqrt_degree[row];

            let mut accumulated = vec![0.0f32; n_components];
            for entry in start..end {
                let column = graph.indices()[entry] as usize;
                let weight = graph.data()[entry] * scale * inverse_sqrt_degree[column];
                for component in 0..n_components {
                    accumulated[component] += weight * basis[[column, component]];
                }
            }
            for component in 0..n_components {
                accumulated[component] += basis[[row, component]];
            }
            accumulated
        })
        .collect::<Vec<_>>();

    let mut result = Array2::zeros((graph.rows(), n_components));
    for (row, values) in rows.into_iter().enumerate() {
        for (component, value) in values.into_iter().enumerate() {
            result[[row, component]] = value;
        }
    }
    result
}

fn deflate(basis: &mut Array2<f32>, trivial: &Array1<f32>) {
    for mut column in basis.columns_mut() {
        let projection = column
            .iter()
            .zip(trivial.iter())
            .map(|(a, b)| a * b)
            .sum::<f32>();
        column.zip_mut_with(trivial, |value, t| *value -= projection * t);
    }
}

/// Modified Gram-Schmidt, sequential over columns so the result does not depend on
/// thread scheduling.
fn orthonormalise(basis: &mut Array2<f32>) {
    for component in 0..basis.ncols() {
        for earlier in 0..component {
            let projection = basis
                .column(component)
                .iter()
                .zip(basis.column(earlier).iter())
                .map(|(a, b)| a * b)
                .sum::<f32>();
            let subtracted = basis.column(earlier).mapv(|value| value * projection);
            basis
                .column_mut(component)
                .zip_mut_with(&subtracted, |value, s| *value -= s);
        }

        let norm = basis
            .column(component)
            .iter()
            .map(|v| v * v)
            .sum::<f32>()
            .sqrt();
        if norm > f32::EPSILON {
            basis
                .column_mut(component)
                .mapv_inplace(|value| value / norm);
        }
    }
}

fn normalise(vector: &mut Array1<f32>) {
    let norm = vector.iter().map(|v| v * v).sum::<f32>().sqrt();
    if norm > f32::EPSILON {
        vector.mapv_inplace(|value| value / norm);
    }
}

/// UMAP expects an embedding roughly spanning [-10, 10]. The jitter breaks ties between
/// points that the eigenvectors place at exactly the same spot.
fn rescale(basis: &mut Array2<f32>, seed: u64) {
    let largest = basis.iter().fold(0.0f32, |acc, value| acc.max(value.abs()));
    if largest > f32::EPSILON {
        basis.mapv_inplace(|value| value * TARGET_SPAN / largest);
    }

    let mut rng = StdRng::seed_from_u64(seed);
    basis.mapv_inplace(|value| value + rng.random_range(-JITTER..JITTER));
}

fn random_basis(n: usize, n_components: usize, seed: u64) -> Array2<f32> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut basis = Array2::from_shape_fn((n, n_components), |_| rng.random_range(-1.0..1.0));
    orthonormalise(&mut basis);
    basis
}
