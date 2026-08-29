//! `spectral_init` takes the top eigenvectors of the graph Laplacian, and a graph of nearly
//! disconnected groups has one near-null direction per group. Asking for fewer dimensions
//! than there are groups leaves the answer undetermined, so the seed picks which projection
//! of a degenerate eigenspace comes back.

use ndarray::Array2;
use rosella::embedding::spectral::spectral_init;
use sprs::{CsMatI, TriMatI};

const PER_BLOCK: usize = 40;
const BLOCKS: usize = 12;

/// Tight blocks joined in a ring, which is the shape a k-NN graph over clustered contigs
/// has: many groups that touch only weakly, so the top eigenvalues sit close together.
fn blocked_graph(blocks: usize) -> CsMatI<f32, u32, usize> {
    let n = blocks * PER_BLOCK;
    let mut triplets = TriMatI::new((n, n));
    let mut add = |i: usize, j: usize, weight: f32| {
        triplets.add_triplet(i, j, weight);
        triplets.add_triplet(j, i, weight);
    };

    for block in 0..blocks {
        let start = block * PER_BLOCK;
        for i in start..start + PER_BLOCK {
            for j in i + 1..start + PER_BLOCK {
                add(i, j, 1.0);
            }
        }
        add(start, (start + PER_BLOCK) % n, 0.01);
    }
    triplets.to_csr()
}

fn unit_columns(mut basis: Array2<f32>) -> Array2<f32> {
    for mut column in basis.columns_mut() {
        let norm = column.iter().map(|v| v * v).sum::<f32>().sqrt();
        if norm > f32::EPSILON {
            column.mapv_inplace(|value| value / norm);
        }
    }
    basis
}

/// Mean squared cosine of the principal angles between two column spaces. One means the two
/// span the same subspace; `n_components / n` is what two random subspaces give.
fn agreement_across_seeds(graph: &CsMatI<f32, u32, usize>, n_components: usize) -> f64 {
    let left = unit_columns(spectral_init(graph, n_components, 42));
    let right = unit_columns(spectral_init(graph, n_components, 7));

    let mut total = 0.0;
    for a in left.columns() {
        for b in right.columns() {
            let dot = a.iter().zip(b.iter()).map(|(x, y)| x * y).sum::<f32>() as f64;
            total += dot * dot;
        }
    }
    total / n_components as f64
}

#[test]
fn the_same_seed_gives_the_same_subspace() {
    let graph = blocked_graph(BLOCKS);
    assert_eq!(spectral_init(&graph, 5, 42), spectral_init(&graph, 5, 42));
}

/// The measured shape, and the reason the embedding dimensionality and the run to run spread
/// are one problem rather than two: too few dimensions and the seed decides which ones.
#[test]
fn too_few_dimensions_leave_the_seed_to_choose_them() {
    let graph = blocked_graph(BLOCKS);
    let under = agreement_across_seeds(&graph, 5);
    let matched = agreement_across_seeds(&graph, BLOCKS);

    assert!(
        matched > under + 0.3,
        "asking for one dimension per group agreed {matched:.4} and asking for five \
         agreed {under:.4}, so the degeneracy is not what decides the subspace"
    );
}
