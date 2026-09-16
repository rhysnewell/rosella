//! Neither arm holds both planted genomes whole, so an assembly-wide choice between them cannot
//! recover both and a per-bin arbitration is the only thing that can.

use std::collections::{HashMap, HashSet};

use rosella::clustering::clusterer::{Partitioning, find_partitions};
use rosella::clustering::graph_partition::Partition;
use rosella::quality::{Quality, Scorer};
use rosella::recover::ladder::{Judge, best_per_arm, combine};

#[path = "../support/bars.rs"]
mod bars;
use sprs::{CsMatI, TriMatI};

const GENOME: usize = 4;
const CONTIGS: usize = 8;

struct Planted;

impl Scorer for Planted {
    fn score(&self, contigs: &[usize]) -> Quality {
        let mut per_genome = [0usize; 2];
        for contig in contigs {
            per_genome[contig / GENOME] += 1;
        }
        let whole = per_genome[0].max(per_genome[1]);
        let foreign = per_genome[0] + per_genome[1] - whole;
        Quality {
            completeness: 100.0 * whole as f64 / GENOME as f64,
            contamination: 100.0 * foreign as f64 / GENOME as f64,
        }
    }

    fn features(&self, contigs: &[usize]) -> HashSet<u32> {
        contigs.iter().map(|contig| *contig as u32).collect()
    }
}

fn partitioning(clusters: &[&[usize]]) -> Partitioning {
    Partitioning {
        cluster_map: clusters
            .iter()
            .enumerate()
            .map(|(id, contigs)| (id, contigs.iter().copied().collect::<HashSet<_>>()))
            .collect::<HashMap<_, _>>(),
        outliers: HashSet::new(),
        score: None,
        arm: Partition::Leiden,
        seed: 0,
    }
}

fn bins(held: &Partitioning) -> Vec<Vec<usize>> {
    let mut found = held
        .cluster_map
        .values()
        .map(|contigs| {
            let mut contigs = contigs.iter().copied().collect::<Vec<_>>();
            contigs.sort_unstable();
            contigs
        })
        .collect::<Vec<_>>();
    found.sort_unstable();
    found
}

fn combined() -> Partitioning {
    let leiden = partitioning(&[&[0, 1, 2, 3], &[4, 5], &[6, 7]]);
    let labelprop = partitioning(&[&[0, 1], &[2, 3], &[4, 5, 6, 7]]);
    let contigs = (0..CONTIGS).collect::<Vec<_>>();
    let scorer = Planted;
    let judge = Judge {
        quality: &scorer,
        contigs: &contigs,
        bars: bars::bars(80.0),
    };
    combine(vec![leiden, labelprop], &judge)
}

#[test]
fn the_combination_takes_a_whole_genome_from_each_arm() {
    assert_eq!(
        bins(&combined()),
        vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7]],
        "each arm holds one planted genome whole and splits the other, so combining should \
         return both and neither arm's split halves"
    );
}

#[test]
fn every_contig_is_placed_once() {
    let held = combined();
    let mut placed = held
        .cluster_map
        .values()
        .flatten()
        .copied()
        .chain(held.outliers.iter().copied())
        .collect::<Vec<_>>();
    let seen = placed.iter().copied().collect::<HashSet<_>>();
    assert_eq!(placed.len(), seen.len(), "a contig was placed twice");
    placed.sort_unstable();
    assert_eq!(placed, (0..CONTIGS).collect::<Vec<_>>());
}

fn blocked_graph() -> CsMatI<f32, u32, usize> {
    let per_block = 30;
    let n = 4 * per_block;
    let mut triplets = TriMatI::new((n, n));
    let mut add = |i: usize, j: usize, weight: f32| {
        triplets.add_triplet(i, j, weight);
        triplets.add_triplet(j, i, weight);
    };
    for block in 0..4 {
        let start = block * per_block;
        for i in start..start + per_block {
            for j in i + 1..start + per_block {
                add(i, j, 1.0);
            }
        }
        add(start, (start + per_block) % n, 0.01);
    }
    triplets.to_csr()
}

#[test]
fn both_arms_reach_the_ladder_the_pool_reads() {
    let graph = blocked_graph();
    let lengths = vec![10_000; graph.rows()];
    let rungs = |kind| {
        find_partitions(&graph, &lengths, None, 42, kind, true)
            .expect("the ladder is never empty")
            .len()
    };
    assert_eq!(rungs(Partition::LabelProp), 1);
    assert_eq!(rungs(Partition::Both), rungs(Partition::Leiden) + 1);
}

/// The fine Leiden rungs are where the fragments that cut a whole genome come from, so the arms
/// source has to keep one labelling per arm and drop the rest of the ladder.
#[test]
fn the_arms_source_keeps_the_objectives_pick_from_each_arm() {
    let scored = |clusters: &[&[usize]], score: f64, arm: Partition| Partitioning {
        score: Some(score),
        arm,
        ..partitioning(clusters)
    };
    let ladder = vec![
        scored(&[&[0, 1]], 0.9, Partition::Leiden),
        scored(&[&[4, 5, 6, 7]], 0.8, Partition::LabelProp),
        scored(&[&[0, 1, 2, 3], &[4, 5, 6, 7]], 0.7, Partition::Leiden),
        scored(&[&[4, 5]], 0.6, Partition::LabelProp),
    ];

    let contigs = (0..CONTIGS).collect::<Vec<_>>();
    let scorer = Planted;
    let judge = Judge {
        quality: &scorer,
        contigs: &contigs,
        bars: bars::bars(80.0),
    };

    let kept = best_per_arm(ladder, &judge);

    assert_eq!(kept.len(), 2);
    assert_eq!(
        kept.iter()
            .map(|held| held.cluster_map.len())
            .collect::<Vec<_>>(),
        vec![2, 1],
        "each arm sends the rung its markers pick, not the one the objective ranks first"
    );
}

fn judged(arms: Vec<Partitioning>) -> Partitioning {
    let contigs = (0..CONTIGS).collect::<Vec<_>>();
    let scorer = Planted;
    let judge = Judge {
        quality: &scorer,
        contigs: &contigs,
        bars: bars::bars(80.0),
    };
    combine(arms, &judge)
}

/// A candidate that loses part of itself is rescored on what is left and has to win again, which
/// is how a fused proposal still yields the clean genome inside it.
#[test]
fn a_part_claimed_candidate_comes_back_on_what_is_left() {
    let arms = vec![
        partitioning(&[&[0, 1, 2, 3, 4]]),
        partitioning(&[&[4, 5, 6, 7]]),
    ];
    let held = judged(arms);
    assert_eq!(
        bins(&held),
        vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7]],
        "the clean genome claims contig 4, and the remnant is a whole genome in its own right"
    );
    assert!(held.outliers.is_empty());
}
