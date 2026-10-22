//! Counting communities at 50 per cent counts a halved genome twice, which is what took the judge
//! to the 875 community rung on cami_i_high where the truth is 1075 genomes and the coarse rung
//! scores 35 more bins.

use std::collections::HashSet;

use rosella::clustering::clusterer::Partitioning;
use rosella::clustering::graph_partition::Partition;
use rosella::quality::{Quality, Scorer};
use rosella::recover::ladder::{Judge, top_rungs};
use rosella::refine::ranking::sorted;

use crate::bars;

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
            ..Default::default()
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
            .map(|(id, contigs)| (id, sorted(contigs.iter().copied())))
            .collect(),
        outliers: Vec::new(),
        score: None,
        arm: Partition::Leiden,
        seed: 0,
    }
}

fn chosen(completeness: f64) -> Vec<usize> {
    let contigs = (0..CONTIGS).collect::<Vec<_>>();
    let scorer = Planted;
    let judge = Judge {
        quality: &scorer,
        contigs: &contigs,
        bars: bars::bars(completeness),
    };
    let ladder = vec![
        partitioning(&[&[0, 1, 2, 3], &[4, 5, 6, 7]]),
        partitioning(&[&[0, 1], &[2, 3], &[4, 5], &[6, 7]]),
    ];
    top_rungs(ladder, &judge)
        .iter()
        .map(|held| held.cluster_map.len())
        .collect()
}

#[test]
fn the_judge_counts_at_the_accept_bar_it_is_given() {
    assert_eq!(
        chosen(80.0),
        [2, 4],
        "an 80 per cent bar refuses the halves"
    );
    assert_eq!(
        chosen(50.0),
        [4, 2],
        "a 50 per cent bar counts each half as a genome"
    );
}

struct Pairs;

impl Scorer for Pairs {
    fn score(&self, contigs: &[usize]) -> Quality {
        let genomes = contigs
            .iter()
            .map(|contig| contig / 2)
            .collect::<HashSet<_>>();
        Quality {
            completeness: 50.0 * contigs.len() as f64 / genomes.len() as f64,
            contamination: 100.0 * (genomes.len() - 1) as f64,
            ..Default::default()
        }
    }

    fn features(&self, contigs: &[usize]) -> HashSet<u32> {
        contigs.iter().map(|contig| *contig as u32).collect()
    }
}

#[test]
fn the_runner_up_rung_is_kept_and_the_third_is_not() {
    let contigs = (0..10).collect::<Vec<_>>();
    let judge = Judge {
        quality: &Pairs,
        contigs: &contigs,
        bars: bars::bars(80.0),
    };
    let ladder = vec![
        partitioning(&[&[0, 1], &[2, 3], &[4, 5], &[6, 7], &[8], &[9]]),
        partitioning(&[&[0, 1], &[2, 3], &[4, 5], &[6], &[7], &[8], &[9]]),
        partitioning(&[&[0, 1], &[2, 3], &[4], &[5], &[6], &[7], &[8], &[9]]),
    ];
    let kept = top_rungs(ladder, &judge)
        .iter()
        .map(|held| held.cluster_map.len())
        .collect::<Vec<_>>();
    assert_eq!(
        kept,
        [6, 7],
        "the rungs with four and three bins over the bar go, the one with two does not"
    );
}
