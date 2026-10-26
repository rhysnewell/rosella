//! Blocks are planted with random edges rather than as cliques, because a clique fixture
//! gives refinement nothing to choose between and the sampler cannot show itself.

use crate::leiden;
use rayon::prelude::*;
use rosella::clustering::graph_partition::label_propagation;
use rosella::clustering::layout::{base_level, base_level_in};
use rosella::clustering::leiden::{leiden_from, resolutions};
use sprs::{CsMatI, TriMatI};

const PER_BLOCK: usize = 25;
const BLOCKS: usize = 8;
const WITHIN: f64 = 0.35;
const ACROSS: f64 = 0.02;

fn planted_graph() -> CsMatI<f32, u32, usize> {
    let n = BLOCKS * PER_BLOCK;
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };

    let mut triplets = TriMatI::new((n, n));
    for i in 0..n {
        for j in i + 1..n {
            let probability = if i / PER_BLOCK == j / PER_BLOCK {
                WITHIN
            } else {
                ACROSS
            };
            if next() < probability {
                let weight = (0.5 + next()) as f32;
                triplets.add_triplet(i, j, weight);
                triplets.add_triplet(j, i, weight);
            }
        }
    }
    triplets.to_csr()
}

#[test]
fn the_ladder_is_independent_of_the_rayon_schedule() {
    let graph = planted_graph();
    let ladder = resolutions(&graph, &vec![1.0; graph.rows()], 8);

    let serial = ladder
        .iter()
        .map(|resolution| leiden(&graph, *resolution, 42))
        .collect::<Vec<_>>();
    let parallel = ladder
        .par_iter()
        .map(|resolution| leiden(&graph, *resolution, 42))
        .collect::<Vec<_>>();

    assert_eq!(serial, parallel);
}

// Unit weights make gains tie and self loops make a node read itself, so a layout that broke ties
// or staleness on its own numbering would show.
#[test]
fn a_base_laid_out_in_any_order_labels_as_the_given_order() {
    let planted = planted_graph();
    let mut tied = planted.clone();
    tied.data_mut().fill(1.0);
    let looped = &tied + &CsMatI::<f32, u32, usize>::eye(tied.rows());
    for graph in [planted, tied, looped] {
        let sizes = vec![1.0; graph.rows()];
        let given = (0..graph.rows() as u32).collect::<Vec<_>>();
        let reversed = given.iter().rev().copied().collect::<Vec<_>>();
        let wanted = base_level_in(&graph, &sizes, given);
        for base in [
            base_level(&graph, &sizes),
            base_level_in(&graph, &sizes, reversed),
        ] {
            assert_eq!(label_propagation(&base, 7), label_propagation(&wanted, 7));
            for resolution in resolutions(&graph, &sizes, 4) {
                let laid = leiden_from(&base, resolution, 7);
                assert_eq!(laid, leiden_from(&wanted, resolution, 7), "{resolution}");
            }
        }
    }
}
