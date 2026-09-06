//! Blocks are planted with random edges rather than as cliques, because a clique fixture
//! gives refinement nothing to choose between and the sampler cannot show itself.

use rayon::prelude::*;
use rosella::clustering::leiden::{leiden, resolutions};
use sprs::{CsMatI, TriMatI};
use std::collections::HashSet;

const PER_BLOCK: usize = 25;
const BLOCKS: usize = 8;
const WITHIN: f64 = 0.35;
const ACROSS: f64 = 0.02;
const THETA: f64 = 0.5;

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

fn planted(labels: &[i32]) -> f64 {
    let mut agree = 0.0;
    let mut total = 0.0;
    for i in 0..labels.len() {
        for j in i + 1..labels.len() {
            if (i / PER_BLOCK == j / PER_BLOCK) == (labels[i] == labels[j]) {
                agree += 1.0;
            }
            total += 1.0;
        }
    }
    agree / total
}

fn communities(labels: &[i32]) -> usize {
    labels.iter().collect::<HashSet<_>>().len()
}

#[test]
fn sampling_is_independent_of_the_rayon_schedule() {
    let graph = planted_graph();
    let ladder = resolutions(&graph, None, 8);

    let serial = ladder
        .iter()
        .map(|resolution| leiden(&graph, None, *resolution, Some(THETA), 42))
        .collect::<Vec<_>>();
    let parallel = ladder
        .par_iter()
        .map(|resolution| leiden(&graph, None, *resolution, Some(THETA), 42))
        .collect::<Vec<_>>();

    assert_eq!(serial, parallel);
}

#[test]
fn sampling_reaches_a_partition_the_argmax_cannot() {
    let graph = planted_graph();

    let differs = resolutions(&graph, None, 8).iter().any(|resolution| {
        let greedy = leiden(&graph, None, *resolution, None, 42);
        (42..48).any(|seed| leiden(&graph, None, *resolution, Some(THETA), seed) != greedy)
    });

    assert!(
        differs,
        "no resolution and seed moved off the greedy labelling, so theta is being ignored"
    );
}

#[test]
fn sampling_still_recovers_the_planted_blocks() {
    let graph = planted_graph();
    let best = resolutions(&graph, None, 8)
        .iter()
        .map(|resolution| {
            let labels = leiden(&graph, None, *resolution, Some(THETA), 42);
            (planted(&labels), communities(&labels))
        })
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .expect("the resolution ladder is never empty");

    assert!(
        best.0 > 0.98,
        "the best of the sampled resolutions agreed on only {:.4} of pairs, over {} \
         communities against {BLOCKS} planted",
        best.0,
        best.1
    );
}
