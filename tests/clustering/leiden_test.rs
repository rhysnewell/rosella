//! Blocks are planted with random edges rather than as cliques, because a clique fixture
//! gives refinement nothing to choose between and the sampler cannot show itself.

use rayon::prelude::*;
use rosella::clustering::leiden::{leiden, resolutions};
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
    let ladder = resolutions(&graph, None, 8);

    let serial = ladder
        .iter()
        .map(|resolution| leiden(&graph, None, *resolution, 42))
        .collect::<Vec<_>>();
    let parallel = ladder
        .par_iter()
        .map(|resolution| leiden(&graph, None, *resolution, 42))
        .collect::<Vec<_>>();

    assert_eq!(serial, parallel);
}

