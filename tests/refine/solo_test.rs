//! Whether a bin holding more than one genome-sized contig comes apart into them, and
//! where genome-sized is read from.

use std::collections::BTreeMap;

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::homology::{Homology, HomologySettings, Pair};
use rosella::refine::solo::{SoloPool, candidate, floor};

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

    assert_eq!(floor(&features, &bins, &[3], MIN_BIN_SIZE, SoloPool::Alone), Some(1_300_000));
    assert_eq!(floor(&features, &BTreeMap::new(), &[2], MIN_BIN_SIZE, SoloPool::Alone), None);
}

/// A closed genome that collected short contigs still carries most of its bin, and a bin of
/// two equal genome-sized contigs carries neither, so only `Long` counts both halves.
#[test]
fn a_genome_holding_most_of_its_bin_measures_the_floor_when_nothing_stands_alone() {
    let lengths = [3_000_000, 10_000, 10_000, 2_400_000, 2_400_000, 500_000, 20_000];
    let (coverage, tnf) = features(&lengths);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let bins = BTreeMap::from([(0, vec![0, 1, 2]), (1, vec![3, 4]), (2, vec![5, 6])]);

    assert_eq!(floor(&features, &bins, &[], MIN_BIN_SIZE, SoloPool::Alone), None);
    assert_eq!(floor(&features, &bins, &[], MIN_BIN_SIZE, SoloPool::Majority), Some(1_500_000));
    assert_eq!(floor(&features, &bins, &[], MIN_BIN_SIZE, SoloPool::Long), Some(1_200_000));
}

#[test]
fn two_genome_sized_contigs_each_stand_alone_and_the_rest_stays_together() {
    let lengths = [3_000_000, 10_000, 2_400_000, 10_000, 10_000];
    let (coverage, tnf) = features(&lengths);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);

    let kept = candidate(&features, &[0, 1, 2, 3, 4], 1_300_000, false).expect("a split");
    assert_eq!(kept, vec![vec![0], vec![2], vec![1, 3, 4]]);
    assert!(candidate(&features, &[0, 1, 3, 4], 1_300_000, false).is_none());
}

/// Length alone breaks one genome in two roughly one firing in five. Homology is the evidence
/// that says which case this is: two loci of one genome do not align to each other.
#[test]
fn homology_decides_whether_genome_sized_contigs_come_apart() {
    let lengths = [3_000_000, 10_000, 2_400_000, 10_000, 10_000];
    let (coverage, tnf) = features(&lengths);
    let settings = HomologySettings::default();

    let silent = Homology::from_pairs([], &lengths, settings);
    let one_genome = ContigFeatures::new(&coverage, &tnf, &lengths).with_homology(Some(&silent));
    assert!(candidate(&one_genome, &[0, 1, 2, 3, 4], 1_300_000, false).is_none());

    let aligned = Homology::from_pairs(
        [Pair {
            one: 0,
            other: 2,
            identity: 98.0,
            aligned_one: 70.0,
            aligned_other: 65.0,
        }],
        &lengths,
        settings,
    );
    let two_organisms =
        ContigFeatures::new(&coverage, &tnf, &lengths).with_homology(Some(&aligned));
    let kept = candidate(&two_organisms, &[0, 1, 2, 3, 4], 1_300_000, false).expect("a split");
    assert_eq!(kept, vec![vec![0], vec![2], vec![1, 3, 4]]);
}

fn a_long_pair_and_fragments() -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let depth = [10.0, 40.0, 10.0, 25.0, 25.0];
    let shape = [
        [4.0, 1.0, 1.0, 1.0],
        [1.0, 4.0, 1.0, 1.0],
        [4.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 4.0, 1.0],
        [1.0, 1.0, 4.0, 1.0],
    ];
    let mut coverage = Array2::zeros((5, 2));
    let mut tnf = Array2::zeros((5, 4));
    for row in 0..5 {
        coverage[[row, 0]] = depth[row];
        coverage[[row, 1]] = 4.0;
        for column in 0..4 {
            tnf[[row, column]] = shape[row][column];
        }
    }
    (
        coverage,
        tnf,
        vec![3_000_000, 2_400_000, 250_000, 250_000, 250_000],
    )
}

/// Contig 2 belongs with the long contig it was assembled beside, and one leftover bin holding
/// every short contig strands it away from its own genome.
#[test]
fn a_short_contig_follows_the_piece_it_sits_nearest() {
    let (coverage, tnf, lengths) = a_long_pair_and_fragments();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);

    let swept = candidate(&features, &[0, 1, 2, 3, 4], 1_300_000, false).expect("a split");
    let scattered = candidate(&features, &[0, 1, 2, 3, 4], 1_300_000, true).expect("a split");

    assert_eq!(swept, vec![vec![0], vec![1], vec![2, 3, 4]]);
    assert_eq!(scattered, vec![vec![0, 2], vec![1], vec![3, 4]]);
}
