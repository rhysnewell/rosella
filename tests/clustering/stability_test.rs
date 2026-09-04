//! Stability is only a selector if it can tell a real partition from a coarse one. These pin
//! the chance correction working and the one case it does not cover.

use rosella::clustering::stability::adjusted_rand_index;

const NODES: usize = 120;
const GROUPS: usize = 6;

fn blocks() -> Vec<i32> {
    (0..NODES)
        .map(|node| (node / (NODES / GROUPS)) as i32)
        .collect()
}

fn interleaved() -> Vec<i32> {
    (0..NODES)
        .map(|node| ((node * 7) % GROUPS) as i32)
        .collect()
}

#[test]
fn agreement_is_corrected_for_chance() {
    let blocks = blocks();
    let self_agreement = adjusted_rand_index(&blocks, &blocks);
    assert!(
        (self_agreement - 1.0).abs() < 1e-12,
        "a labelling should agree with itself, got {self_agreement}"
    );

    let scrambled = adjusted_rand_index(&blocks, &interleaved());
    assert!(
        scrambled.abs() < 0.05,
        "an unrelated labelling should score near zero, got {scrambled}"
    );
}

#[test]
fn two_coarse_partitions_agree_for_free() {
    let one = vec![0i32; NODES];
    let agreement = adjusted_rand_index(&one, &one);
    assert!(
        (agreement - 1.0).abs() < 1e-12,
        "the trivial partition is perfectly stable, got {agreement}, which is why stability \
         needs a null before it can rank a ladder"
    );
}
