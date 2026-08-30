use std::cmp::Ordering;

use ndarray::{Array1, Array2};
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::embedding::Graph;

const POWER_ITERATIONS: usize = 200;
const TARGET_SPAN: f32 = 10.0;
const JITTER: f32 = 1e-4;
const GOLDEN_FRACTION: f64 = 0.618_033_988_749_895;

pub const SPECTRAL_INIT_NAMES: [&str; 2] = ["random", "landmark"];

/// Which subspace power iteration converges to, when the graph has more near-disconnected
/// groups than dimensions asked for and every choice is an equally valid answer.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum SpectralInit {
    /// A uniform random start, which leaves the choice to the seed.
    #[default]
    Random,
    /// Diffusion from graph landmarks. Seed free, and measured at a third of the anchor's
    /// `t1` on CAMI I high: ten localised bumps separate ten regions and mix none of the rest.
    Landmark,
}

impl SpectralInit {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "landmark" => Some(Self::Landmark),
            "random" => Some(Self::Random),
            _ => None,
        }
    }
}

/// Eigenvectors of the normalised graph Laplacian, which is what UMAP initialises from.
/// A uniform random start leaves flight's low `b` clumping the embedding into hundreds of
/// disconnected knots that the optimiser never merges.
pub fn spectral_init(
    graph: &Graph,
    n_components: usize,
    seed: u64,
    init: SpectralInit,
) -> Array2<f32> {
    let n = graph.rows();
    if n == 0 {
        return Array2::zeros((0, n_components));
    }

    let degrees = degrees(graph);
    let inverse_sqrt_degree = degrees.mapv(|d| if d > 0.0 { 1.0 / d.sqrt() } else { 0.0 });

    // The all-ones vector scaled by the degrees is the Laplacian's trivial eigenvector.
    // It carries no layout information, so it gets projected out every iteration.
    let mut trivial = degrees.mapv(|d| d.sqrt());
    normalise(&mut trivial);

    let mut basis = match init {
        SpectralInit::Landmark => {
            landmark_basis(graph, &inverse_sqrt_degree, &degrees, n_components)
        }
        SpectralInit::Random => random_basis(n, n_components, seed),
    };
    for _ in 0..POWER_ITERATIONS {
        basis = shifted_multiply(graph, &inverse_sqrt_degree, &basis);
        deflate(&mut basis, &trivial);
        orthonormalise(&mut basis);
    }

    rescale(&mut basis, seed, init);
    basis
}

/// Columns sharing a quotient span a degenerate eigenspace, which is the state that leaves
/// the subspace undetermined and the seed free to pick one. Logged so that state is visible
/// on a real assembly rather than only on the synthetic graph the tests use.
pub fn rayleigh_quotients(graph: &Graph, basis: &Array2<f32>) -> Vec<f32> {
    let inverse_sqrt_degree = degrees(graph).mapv(|d| if d > 0.0 { 1.0 / d.sqrt() } else { 0.0 });
    let product = shifted_multiply(graph, &inverse_sqrt_degree, basis);
    let mut quotients = (0..basis.ncols())
        .map(|column| {
            basis
                .column(column)
                .iter()
                .zip(product.column(column).iter())
                .map(|(a, b)| a * b)
                .sum::<f32>()
        })
        .collect::<Vec<_>>();
    quotients.sort_by(|a, b| b.total_cmp(a));
    quotients
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
fn rescale(basis: &mut Array2<f32>, seed: u64, init: SpectralInit) {
    let largest = basis.iter().fold(0.0f32, |acc, value| acc.max(value.abs()));
    if largest > f32::EPSILON {
        basis.mapv_inplace(|value| value * TARGET_SPAN / largest);
    }

    match init {
        SpectralInit::Landmark => {
            let columns = basis.ncols();
            for ((row, column), value) in basis.indexed_iter_mut() {
                *value += golden_offset(row * columns + column);
            }
        }
        SpectralInit::Random => {
            let mut rng = StdRng::seed_from_u64(seed);
            basis.mapv_inplace(|value| value + rng.random_range(-JITTER..JITTER));
        }
    }
}

/// Successive multiples of the golden ratio fall as far from each other as an equidistributed
/// sequence can, so distinct positions separate without a seed.
fn golden_offset(position: usize) -> f32 {
    let fraction = ((position + 1) as f64 * GOLDEN_FRACTION).fract();
    JITTER * (2.0 * fraction as f32 - 1.0)
}

/// Diffusion from landmarks chosen for being least reachable from those already taken. The
/// subspace this converges to is the data's own, where a random start leaves the seed to pick
/// one of the many the degenerate eigenspace admits.
fn landmark_basis(
    graph: &Graph,
    inverse_sqrt_degree: &Array1<f32>,
    degrees: &Array1<f32>,
    n_components: usize,
) -> Array2<f32> {
    let n = graph.rows();
    let steps = diffusion_steps(graph);
    let mut reach = Array1::<f32>::zeros(n);
    let mut basis = Array2::<f32>::zeros((n, n_components));

    for component in 0..n_components {
        let landmark = next_landmark(&reach, degrees);
        let mut mass = Array2::<f32>::zeros((n, 1));
        mass[[landmark, 0]] = 1.0;

        for _ in 0..steps {
            mass = shifted_multiply(graph, inverse_sqrt_degree, &mass);
            let norm = mass.iter().map(|v| v * v).sum::<f32>().sqrt();
            if norm > f32::EPSILON {
                mass.mapv_inplace(|value| value / norm);
            }
        }

        for row in 0..n {
            basis[[row, component]] = mass[[row, 0]];
            reach[row] += mass[[row, 0]].abs();
        }
    }

    orthonormalise(&mut basis);
    basis
}

/// Hops for mass to cross the graph, which is what it takes for `reach` to say anything about
/// the far side of it. Derived from the graph rather than pinned, because a denser graph needs
/// fewer and the k-NN degree is set by `--n-neighbours`.
fn diffusion_steps(graph: &Graph) -> usize {
    let n = graph.rows();
    if n < 2 {
        return 1;
    }
    let mean_degree = graph.nnz() as f64 / n as f64;
    if mean_degree <= 1.0 {
        return 1;
    }
    ((n as f64).ln() / mean_degree.ln()).ceil() as usize
}

/// Least reached, then most connected, then lowest index. The degree tie-break matters where
/// nothing has been reached yet and where a component is isolated: it takes the hub of an
/// untouched region rather than a stray singleton.
fn next_landmark(reach: &Array1<f32>, degrees: &Array1<f32>) -> usize {
    let mut best = 0usize;
    for row in 1..reach.len() {
        let closer = reach[row].total_cmp(&reach[best]);
        if closer == Ordering::Less
            || (closer == Ordering::Equal && degrees[row] > degrees[best])
        {
            best = row;
        }
    }
    best
}

fn random_basis(n: usize, n_components: usize, seed: u64) -> Array2<f32> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut basis = Array2::from_shape_fn((n, n_components), |_| rng.random_range(-1.0..1.0));
    orthonormalise(&mut basis);
    basis
}
