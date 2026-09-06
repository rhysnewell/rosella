//! Whether a bin holding more than one genome-sized contig comes apart into them, and
//! where genome-sized is read from.

use std::collections::BTreeMap;

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::solo::{candidate, floor};

const MIN_BIN_SIZE: usize = 200_000;

fn features(lengths: &[usize]) -> (Array2<f64>, Array2<f64>) {
    (
        Array2::zeros((lengths.len(), 2)),
        Array2::zeros((lengths.len(), 6)),
    )
}

/// A lone short contig is an outlier, not a genome, so it says nothing about genome size.
#[test]
fn the_floor_is_half_the_median_closed_genome() {
    let lengths = [3_000_000, 2_000_000, 50_000, 2_600_000, 20_000, 10_000];
    let (coverage, tnf) = features(&lengths);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let bins = BTreeMap::from([(0, vec![0]), (1, vec![1]), (2, vec![2]), (3, vec![4, 5])]);

    assert_eq!(floor(&features, &bins, &[3], MIN_BIN_SIZE), Some(1_300_000));
    assert_eq!(floor(&features, &BTreeMap::new(), &[2], MIN_BIN_SIZE), None);
}

#[test]
fn two_genome_sized_contigs_each_stand_alone_and_the_rest_stays_together() {
    let lengths = [3_000_000, 10_000, 2_400_000, 10_000, 10_000];
    let (coverage, tnf) = features(&lengths);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);

    let kept = candidate(&features, &[0, 1, 2, 3, 4], 1_300_000).expect("a split");
    assert_eq!(kept, vec![vec![0], vec![2], vec![1, 3, 4]]);
    assert!(candidate(&features, &[0, 1, 3, 4], 1_300_000).is_none());
}
