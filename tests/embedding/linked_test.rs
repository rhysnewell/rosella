//! A link is a coin flip as a must-link on the one real set with a gold, so it may only raise an
//! edge, never replace the neighbour graph, and it has to survive the subset renumbering.

use rosella::embedding::{Graph, linked};
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
fn a_link_raises_a_weak_edge_and_leaves_a_strong_one_alone() {
    let whole = graph(3, &[(0, 1, 0.1), (1, 2, 0.9)]);

    let dense = linked(whole, &[(0, 1), (1, 2)], &[0, 1, 2], 0.5).to_dense();

    assert_eq!(dense[[0, 1]], 0.5);
    assert_eq!(dense[[1, 2]], 0.9);
}

#[test]
fn a_link_between_contigs_with_no_edge_becomes_one() {
    let whole = graph(3, &[(0, 1, 0.1)]);

    let dense = linked(whole, &[(0, 2)], &[0, 1, 2], 0.4).to_dense();

    assert_eq!(dense[[0, 2]], 0.4);
    assert_eq!(dense[[2, 0]], 0.4);
}

#[test]
fn links_are_read_in_the_subsets_own_numbering() {
    let whole = graph(2, &[(0, 1, 0.1)]);

    let dense = linked(whole, &[(3, 7), (3, 4)], &[3, 7], 0.6).to_dense();

    assert_eq!(
        dense[[0, 1]],
        0.6,
        "the pair inside the subset renumbers to 0 and 1"
    );
    assert_eq!(dense.iter().filter(|weight| **weight > 0.0).count(), 2);
}
