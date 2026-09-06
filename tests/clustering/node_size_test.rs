//! Two closed genomes sit beside two fragmented ones. Counting contigs, no rung of the ladder
//! keeps every genome whole and clean in bases at once. Counting bases, several do.

use rosella::clustering::graph_partition::NodeSize;
use rosella::clustering::leiden::{leiden, resolutions};
use sprs::{CsMatI, TriMatI};
use std::collections::HashMap;

const FRAGMENTS: usize = 30;
const FRAGMENT_BP: usize = 2_000;
const CLOSED_BP: usize = 3_000_000;
const WITHIN: f64 = 0.35;
const ACROSS: f64 = 0.02;
const FOREIGN_EDGES: usize = 6;
const BP_TIER: f64 = 0.9;

fn planted() -> (CsMatI<f32, u32, usize>, Vec<usize>, Vec<usize>) {
    let n = 2 * FRAGMENTS + 2;
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };

    let genome = (0..n)
        .map(|i| {
            if i < 2 * FRAGMENTS {
                i / FRAGMENTS
            } else {
                i - 2 * FRAGMENTS + 2
            }
        })
        .collect::<Vec<_>>();
    let lengths = (0..n)
        .map(|i| {
            if i < 2 * FRAGMENTS {
                FRAGMENT_BP
            } else {
                CLOSED_BP
            }
        })
        .collect::<Vec<_>>();

    let mut triplets = TriMatI::new((n, n));
    let mut link = |a: usize, b: usize, weight: f32| {
        triplets.add_triplet(a, b, weight);
        triplets.add_triplet(b, a, weight);
    };
    for i in 0..2 * FRAGMENTS {
        for j in i + 1..2 * FRAGMENTS {
            let probability = if genome[i] == genome[j] {
                WITHIN
            } else {
                ACROSS
            };
            if next() < probability {
                link(i, j, (0.5 + next()) as f32);
            }
        }
    }
    let closed = [2 * FRAGMENTS, 2 * FRAGMENTS + 1];
    link(closed[0], closed[1], 1.0);
    for (slot, node) in closed.iter().enumerate() {
        for k in 0..FOREIGN_EDGES {
            link(*node, slot * FRAGMENTS + k, 1.0);
        }
    }
    (triplets.to_csr(), lengths, genome)
}

fn every_genome_clean(labels: &[i32], lengths: &[usize], genome: &[usize]) -> bool {
    let mut community_bp: HashMap<i32, usize> = HashMap::new();
    let mut shared: HashMap<(usize, i32), usize> = HashMap::new();
    let mut genome_bp: HashMap<usize, usize> = HashMap::new();
    for (node, label) in labels.iter().enumerate() {
        *community_bp.entry(*label).or_default() += lengths[node];
        *shared.entry((genome[node], *label)).or_default() += lengths[node];
        *genome_bp.entry(genome[node]).or_default() += lengths[node];
    }
    genome_bp.iter().all(|(g, total)| {
        let (label, bp) = shared
            .iter()
            .filter(|((owner, _), _)| owner == g)
            .map(|((_, label), bp)| (*label, *bp))
            .max_by_key(|(_, bp)| *bp)
            .expect("every genome owns a node");
        bp as f64 >= BP_TIER * *total as f64 && bp as f64 >= BP_TIER * community_bp[&label] as f64
    })
}

fn rungs_that_get_every_genome(node_size: NodeSize) -> usize {
    let (graph, lengths, genome) = planted();
    let sized = node_size.apply(&graph, &lengths);
    let sizes = sized.sizes.as_deref();
    resolutions(&sized.graph, sizes, 10)
        .iter()
        .filter(|gamma| {
            every_genome_clean(
                &leiden(&sized.graph, sizes, **gamma, None, 42),
                &lengths,
                &genome,
            )
        })
        .count()
}

#[test]
fn bases_separate_closed_genomes_where_contig_counts_cannot() {
    assert_eq!(rungs_that_get_every_genome(NodeSize::Count), 0);
    assert!(rungs_that_get_every_genome(NodeSize::Bp) > 0);
}

#[test]
fn count_sizing_leaves_the_ladder_unchanged() {
    let (graph, lengths, _) = planted();
    let sized = NodeSize::Count.apply(&graph, &lengths);
    assert!(sized.sizes.is_none());
    assert_eq!(
        resolutions(&sized.graph, None, 6),
        resolutions(&graph, None, 6)
    );
}
