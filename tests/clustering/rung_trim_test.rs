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
fn only_a_rung_that_keeps_every_marker_contig_together_trims_a_bin() {
    let splits_the_markers = partitioning(&[&[0, 1, 3, 4], &[2], &[5]]);
    let keeps_the_markers = partitioning(&[&[0, 1, 2, 3], &[4, 5]]);
    let rungs = Rungs::new(&[splits_the_markers, keeps_the_markers], 6);
    let contigs = (0..6).collect::<Vec<_>>();
    let scorer = MarkersOnFirst;
    let judge = Judge {
        quality: &scorer,
        contigs: &contigs,
        lengths: &[1; 6],
        bars: bars::bars(80.0),
        depth: None,
    };
    let trimmed = rungs.trimmed(partitioning(&[&[0, 1, 2, 3, 4], &[5]]), &judge);
    assert_eq!(
        trimmed.cluster_map.values().cloned().collect::<Vec<_>>(),
        vec![vec![0, 1, 2, 3], vec![5]],
        "the rung that splits contig 2 from the other markers says nothing, and the rung that \
         keeps all three together and leaves contig 4 out sends 4 alone to the outliers"
    );
    assert_eq!(trimmed.outliers, vec![4]);
}
