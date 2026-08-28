//! The layout step exists because umap-rs draws negative samples from an OS seeded RNG and
//! writes the embedding from many threads without synchronisation. These pin the property
//! that swap was made for.

use ndarray::Array2;
use rosella::embedding::{
    layout::{LayoutSettings, optimise},
    umap::CurveParams,
};
use sprs::{CsMatI, TriMatI};

const N: usize = 60;
const DIM: usize = 3;

/// Two well separated rings, so the graph has structure the optimiser can pull apart.
fn graph() -> CsMatI<f32, u32, usize> {
    let mut triplets = TriMatI::new((N, N));
    for i in 0..N {
        for step in 1..=3 {
            let j = (i + step) % N;
            let weight = 1.0 / step as f32;
            triplets.add_triplet(i, j, weight);
            triplets.add_triplet(j, i, weight);
        }
    }
    triplets.to_csr()
}

fn init() -> Array2<f32> {
    Array2::from_shape_fn((N, DIM), |(i, d)| {
        ((i * 7 + d * 13) % 41) as f32 / 41.0 * 20.0 - 10.0
    })
}

fn settings(seed: u64) -> LayoutSettings {
    LayoutSettings {
        curve: CurveParams { a: 1.5, b: 0.3 },
        n_epochs: 50,
        seed,
    }
}

#[test]
fn the_same_seed_gives_the_same_embedding() {
    let first = optimise(&graph(), init(), &settings(42), &[]);
    let second = optimise(&graph(), init(), &settings(42), &[]);
    assert_eq!(first, second);
}

#[test]
fn a_different_seed_gives_a_different_embedding() {
    let first = optimise(&graph(), init(), &settings(42), &[]);
    let second = optimise(&graph(), init(), &settings(7), &[]);
    assert_ne!(first, second);
}

#[test]
fn the_embedding_moves_and_stays_finite() {
    let start = init();
    let end = optimise(&graph(), start.clone(), &settings(42), &[]);
    assert!(end.iter().all(|value| value.is_finite()));
    assert_ne!(end, start);
}

#[test]
fn an_empty_graph_leaves_the_initialisation_alone() {
    let empty = TriMatI::<f32, u32>::new((N, N)).to_csr();
    let start = init();
    assert_eq!(optimise(&empty, start.clone(), &settings(42), &[]), start);
}

/// Unit weights have to leave the layout exactly where the unweighted path put it,
/// otherwise the default is no longer the behaviour every recorded baseline was measured on.
#[test]
fn unit_vertex_weights_change_nothing() {
    let uniform = vec![1.0f32; N];
    assert_eq!(
        optimise(&graph(), init(), &settings(42), &uniform),
        optimise(&graph(), init(), &settings(42), &[])
    );
}

#[test]
fn uneven_vertex_weights_move_the_layout() {
    let weights = (0..N)
        .map(|i| if i < N / 2 { 0.25 } else { 4.0 })
        .collect::<Vec<f32>>();
    assert_ne!(
        optimise(&graph(), init(), &settings(42), &weights),
        optimise(&graph(), init(), &settings(42), &[])
    );
}
