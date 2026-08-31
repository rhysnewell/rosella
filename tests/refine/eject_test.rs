//! The one mechanism that takes a contig back out of a bin.

use std::collections::BTreeMap;

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::bin_stats::LevelSource;
use rosella::refine::eject::eject;

const TIGHT: usize = 10;
const STRAY: usize = TIGHT;
const CONTIG_LENGTH: usize = 100_000;
const QUANTILE: f64 = 0.75;

/// One bin of co-abundant contigs of near identical composition, plus a contig that matches
/// neither.
fn fixture() -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let n = TIGHT + 1;
    let mut coverage = Array2::zeros((n, 4));
    let mut tnf = Array2::zeros((n, 4));

    for row in 0..TIGHT {
        let jitter = row as f64 * 0.01;
        coverage[[row, 0]] = 10.0 + jitter;
        coverage[[row, 1]] = 2.0;
        coverage[[row, 2]] = 4.0 + jitter;
        coverage[[row, 3]] = 1.0;
        for column in 0..4 {
            tnf[[row, column]] = [0.1, -0.2, 0.3, -0.4][column] + jitter;
        }
    }

    coverage[[STRAY, 0]] = 500.0;
    coverage[[STRAY, 1]] = 50.0;
    coverage[[STRAY, 2]] = 0.2;
    coverage[[STRAY, 3]] = 0.1;
    for column in 0..4 {
        tnf[[STRAY, column]] = [-2.0, 2.0, -1.5, 1.8][column];
    }

    (coverage, tnf, vec![CONTIG_LENGTH; n])
}

fn run(min_bin_size: usize) -> (Vec<usize>, BTreeMap<usize, Vec<usize>>) {
    let (coverage, tnf, lengths) = fixture();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let mut bins = BTreeMap::from([(0usize, (0..=TIGHT).collect::<Vec<_>>())]);
    let ejected = eject(
        &features,
        &mut bins,
        LevelSource::Flight,
        QUANTILE,
        1.0,
        min_bin_size,
        42,
    );
    (ejected, bins)
}

#[test]
fn the_contig_that_matches_neither_leaves() {
    let (ejected, bins) = run(200_000);

    assert_eq!(ejected, vec![STRAY]);
    assert_eq!(bins[&0], (0..TIGHT).collect::<Vec<_>>());
}

#[test]
fn a_bin_with_no_room_keeps_everything() {
    let (ejected, bins) = run((TIGHT + 1) * CONTIG_LENGTH);

    assert!(ejected.is_empty());
    assert_eq!(bins[&0].len(), TIGHT + 1);
}
