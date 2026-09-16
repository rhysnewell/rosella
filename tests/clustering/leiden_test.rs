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
    let ladder = resolutions(&graph, None, 8, None);

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

/// The anchored ladder has to aim at the same masses on any assembly, and a band wider than
/// the mass on hand has to come back to it rather than ask for a community bigger than the set.
#[test]
fn an_anchored_ladder_aims_at_the_band() {
    let graph = planted_graph();
    let sizes = vec![1_000_000.0; graph.rows()];
    let total = sizes.iter().sum::<f64>();

    let rungs = resolutions(&graph, Some(&sizes), 6, Some((200_000.0, 15_000_000.0)));
    let aimed = rungs.iter().map(|rung| rung / rungs[0]).collect::<Vec<_>>();

    assert!(total > 15_000_000.0);
    assert!((aimed[0] - 1.0).abs() < 1e-9);
    assert!((aimed[5] - 15_000_000.0 / 200_000.0).abs() < 1e-6);

    let clamped = resolutions(&graph, Some(&sizes), 6, Some((200_000.0, total * 4.0)));
    assert_eq!(
        clamped,
        resolutions(&graph, Some(&sizes), 6, Some((200_000.0, total)))
    );
    assert!(clamped[0] < rungs[0]);
}
