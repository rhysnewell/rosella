//! Worth charges contamination against completeness on one scale, so a nearly whole bin carrying
//! a few per cent outranks a clean piece the pool would adopt a rung sooner. Banding the
//! candidates on the pool's ladder first is what puts them back in the order the pool reads.

use std::collections::{HashMap, HashSet};

use rosella::clustering::clusterer::Partitioning;
use rosella::clustering::graph_partition::Partition;
use rosella::quality::{Quality, Scorer};
use rosella::recover::ladder::{Judge, combine};
use rosella::refine::rung::RUNGS;

#[path = "../support/bars.rs"]
mod bars;

struct Table;

impl Scorer for Table {
    fn score(&self, contigs: &[usize]) -> Quality {
        let (completeness, contamination) = match contigs {
            [0, 1, 2, 3] => (95.0, 8.0),
            [0, 1] => (72.0, 0.0),
            _ => (10.0, 0.0),
        };
        Quality {
            completeness,
            contamination,
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

fn bins(rungs: usize) -> Vec<Vec<usize>> {
    let contigs = (0..4).collect::<Vec<_>>();
    let scorer = Table;
    let judge = Judge {
        quality: &scorer,
        contigs: &contigs,
        bars: bars::bars(80.0),
        rungs,
        size_tie: false,
    };
    let ladder = vec![partitioning(&[&[0, 1, 2, 3]]), partitioning(&[&[0, 1]])];
    let mut found = combine(ladder, &judge)
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

#[test]
fn the_band_outranks_worth() {
    assert_eq!(
        bins(0),
        vec![vec![0, 1, 2, 3]],
        "on worth alone the 95 at 8 contamination takes the contigs"
    );
    assert_eq!(
        bins(RUNGS),
        vec![vec![0, 1], vec![2, 3]],
        "the 72 at no contamination passes a rung the 95 at 8 does not"
    );
}

#[test]
fn the_rung_count_decides_which_bands_exist() {
    assert_eq!(
        bins(1),
        vec![vec![0, 1, 2, 3]],
        "neither candidate clears the top rung, so one band holds both and worth decides again"
    );
}
