use std::collections::HashSet;

use rosella::clustering::clusterer::Partitioning;
use rosella::clustering::graph_partition::Partition;
use rosella::quality::{Quality, Scorer};
use rosella::recover::ladder::{Judge, Rungs};

use crate::bars;

const MARKERS: usize = 3;

struct MarkersOnFirst;

impl Scorer for MarkersOnFirst {
    fn score(&self, _contigs: &[usize]) -> Quality {
        Quality::default()
    }

    fn features(&self, contigs: &[usize]) -> HashSet<u32> {
        contigs
            .iter()
            .filter(|contig| **contig < MARKERS)
            .map(|contig| *contig as u32)
            .collect()
    }
}

fn partitioning(clusters: &[&[usize]]) -> Partitioning {
    Partitioning {
        cluster_map: clusters
            .iter()
            .enumerate()
            .map(|(id, contigs)| (id, contigs.to_vec()))
            .collect(),
        outliers: Vec::new(),
        score: None,
        arm: Partition::Leiden,
        seed: 0,
    }
}

#[test]
fn only_a_group_wholly_inside_the_bin_and_apart_from_its_markers_leaves() {
    let splits_the_markers = partitioning(&[&[0, 1, 3, 4], &[2], &[5]]);
    let keeps_too_few_bases = partitioning(&[&[0, 1, 2], &[3, 4], &[5]]);
    let joins_another_bin = partitioning(&[&[0, 1, 2, 3], &[4, 5]]);
    let stands_alone = partitioning(&[&[0, 1, 2, 3], &[4], &[5]]);
    let contigs = (0..6).collect::<Vec<_>>();
    let scorer = MarkersOnFirst;
    let judge = Judge {
        quality: &scorer,
        contigs: &contigs,
        lengths: &[1, 1, 1, 10, 1, 1],
        bars: bars::bars(80.0),
        depth: None,
    };
    let bins = || partitioning(&[&[0, 1, 2, 3, 4], &[5]]);
    let silent = Rungs::new(
        &[splits_the_markers, keeps_too_few_bases, joins_another_bin],
        6,
    );
    assert_eq!(
        silent.trimmed(bins(), &judge).outliers,
        Vec::<usize>::new(),
        "a rung that splits the markers, keeps them with under half the bases, or puts contig 4 \
         with another bin's contig says nothing about contig 4"
    );
    let trimmed = Rungs::new(&[stands_alone], 6).trimmed(bins(), &judge);
    assert_eq!(
        trimmed.cluster_map.values().cloned().collect::<Vec<_>>(),
        vec![vec![0, 1, 2, 3], vec![5]]
    );
    assert_eq!(trimmed.outliers, vec![4]);
}
