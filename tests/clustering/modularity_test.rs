//! DBCV's failure on a graph partition was a shape failure: it rose monotonically as the
//! partition got finer, so it had no interior optimum to find. These pin the shape rather
//! than the value, because the shape is what a ladder search reads.

use rosella::clustering::modularity::modularity;
use sprs::{CsMatI, TriMatI};

const PER_BLOCK: usize = 40;
const BLOCKS: usize = 12;
const BRIDGE: f32 = 0.01;

fn blocked_graph(self_loop: f32) -> CsMatI<f32, u32, usize> {
    let n = BLOCKS * PER_BLOCK;
    let mut triplets = TriMatI::new((n, n));
    let mut add = |i: usize, j: usize, weight: f32| {
        triplets.add_triplet(i, j, weight);
        triplets.add_triplet(j, i, weight);
    };
    for block in 0..BLOCKS {
        let start = block * PER_BLOCK;
        for i in start..start + PER_BLOCK {
            for j in i + 1..start + PER_BLOCK {
                add(i, j, 1.0);
            }
        }
        add(start, (start + PER_BLOCK) % n, BRIDGE);
    }
    if self_loop > 0.0 {
        for i in 0..n {
            triplets.add_triplet(i, i, self_loop);
        }
    }
    triplets.to_csr()
}

fn every_nth_block(divisor: usize) -> Vec<i32> {
    (0..BLOCKS * PER_BLOCK)
        .map(|node| (node / (PER_BLOCK * divisor)) as i32)
        .collect()
}

fn split_blocks(parts: usize) -> Vec<i32> {
    (0..BLOCKS * PER_BLOCK)
        .map(|node| (node / (PER_BLOCK / parts)) as i32)
        .collect()
}

#[test]
fn modularity_peaks_at_the_planted_blocks() {
    let graph = blocked_graph(0.0);
    let planted = modularity(&graph, &every_nth_block(1), 1.0);
    let merged = modularity(&graph, &every_nth_block(2), 1.0);
    let split = modularity(&graph, &split_blocks(2), 1.0);

    assert!(
        planted > merged,
        "planted {planted} should beat merged {merged}"
    );
    assert!(
        planted > split,
        "planted {planted} should beat split {split}"
    );
}

#[test]
fn a_single_community_scores_zero() {
    let graph = blocked_graph(0.0);
    let one = modularity(&graph, &vec![0i32; BLOCKS * PER_BLOCK], 1.0);
    assert!(one.abs() < 1e-9, "one community should score 0, got {one}");
}

#[test]
fn self_loops_do_not_move_the_score() {
    let plain = modularity(&blocked_graph(0.0), &every_nth_block(1), 1.0);
    let looped = modularity(&blocked_graph(5.0), &every_nth_block(1), 1.0);
    assert!(
        (plain - looped).abs() < 1e-9,
        "self loops changed the score, {plain} against {looped}"
    );
}

#[test]
fn resolution_moves_the_peak_finer() {
    // The planted and split partitions cross at gamma 12.3 on this fixture.
    let graph = blocked_graph(0.0);
    let planted = every_nth_block(1);
    let split = split_blocks(2);
    assert!(modularity(&graph, &planted, 1.0) > modularity(&graph, &split, 1.0));
    assert!(modularity(&graph, &split, 20.0) > modularity(&graph, &planted, 20.0));
}
