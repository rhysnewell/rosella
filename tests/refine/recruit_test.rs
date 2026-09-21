//! Taking a cut genome's own contigs back off the bins that hold them.

use std::collections::BTreeMap;

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::embedding::knn::KnnGraph;
use rosella::refine::recruit::{RecruitSettings, recruit};

use crate::scorer::GenomeScorer;

const LENGTH: usize = 100_000;

fn settings() -> RecruitSettings {
    RecruitSettings {
        floor: 70.0,
        confidence: 0.6,
        completeness: 90.0,
        contamination: 5.0,
        max_bin_size: 15_000_000,
        passes: 4,
    }
}

fn depths(means: &[f64]) -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let coverage =
        Array2::from_shape_fn(
            (means.len(), 2),
            |(row, column)| {
                if column == 0 { means[row] } else { 1.0 }
            },
        );
    let tnf = Array2::from_shape_fn((means.len(), 4), |(_, column)| column as f64);
    (coverage, tnf, vec![LENGTH; means.len()])
}

fn complete(count: usize) -> KnnGraph {
    KnnGraph {
        indices: Array2::from_shape_fn((count, count), |(_, column)| column as u32),
        dists: Array2::zeros((count, count)),
    }
}

#[test]
fn takes_back_a_contig_its_own_bin_explains_better() {
    let (coverage, tnf, lengths) = depths(&[
        10.0, 11.0, 12.0, 13.0, 14.0, 12.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0,
    ]);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let genomes = vec![0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1];
    let mut bins = BTreeMap::from([
        (0usize, vec![0, 1, 2, 3, 4]),
        (1usize, vec![5, 6, 7, 8, 9, 10, 11]),
    ]);

    let ledger = recruit(
        &features,
        &GenomeScorer::new(genomes),
        &complete(lengths.len()),
        &mut bins,
        settings(),
    );

    assert_eq!(ledger.taken, 1);
    assert_eq!(bins[&0], vec![0, 1, 2, 3, 4, 5]);
    assert_eq!(bins[&1], vec![6, 7, 8, 9, 10, 11]);
}

#[test]
fn leaves_a_contig_outside_the_bin_it_is_nearest_to() {
    let (coverage, tnf, lengths) = depths(&[
        10.0, 11.0, 12.0, 13.0, 14.0, 30.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0,
    ]);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let genomes = vec![0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1];
    let mut bins = BTreeMap::from([
        (0usize, vec![0, 1, 2, 3, 4]),
        (1usize, vec![5, 6, 7, 8, 9, 10, 11]),
    ]);

    let ledger = recruit(
        &features,
        &GenomeScorer::new(genomes),
        &complete(lengths.len()),
        &mut bins,
        settings(),
    );

    assert_eq!(ledger.taken, 0);
}

#[test]
fn refuses_a_contig_that_contaminates_the_bin_that_wants_it() {
    let (coverage, tnf, lengths) = depths(&[
        10.0, 11.0, 12.0, 13.0, 14.0, 12.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0,
    ]);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let genomes = vec![0, 0, 0, 0, 0, 2, 1, 1, 1, 1, 1, 1];
    let mut bins = BTreeMap::from([
        (0usize, vec![0, 1, 2, 3, 4]),
        (1usize, vec![5, 6, 7, 8, 9, 10, 11]),
    ]);

    let ledger = recruit(
        &features,
        &GenomeScorer::new(genomes),
        &complete(lengths.len()),
        &mut bins,
        settings(),
    );

    assert_eq!(ledger.taken, 0);
}

#[test]
fn refuses_a_contig_whose_markers_the_bin_already_has() {
    let (coverage, tnf, lengths) = depths(&[
        10.0, 11.0, 12.0, 13.0, 14.0, 12.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0,
    ]);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let genomes = vec![0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1];
    let duplicate = vec![7, 7, 7, 7, 7, 7, 9, 9, 9, 9, 9, 9];
    let mut bins = BTreeMap::from([
        (0usize, vec![0, 1, 2, 3, 4]),
        (1usize, vec![5, 6, 7, 8, 9, 10, 11]),
    ]);

    let ledger = recruit(
        &features,
        &GenomeScorer::new(genomes).with_families(duplicate),
        &complete(lengths.len()),
        &mut bins,
        settings(),
    );

    assert_eq!(ledger.taken, 0);
}
