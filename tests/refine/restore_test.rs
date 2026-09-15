//! Putting a dissolved bin back when the pool broke it into pieces that all miss the bar.

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::restore::restore;
use rosella::refine::rung::Rung;

#[path = "../support/scorer.rs"]
mod scorer;

use scorer::BasesScorer;

const PIECE: usize = 100_000;
const CONTIGS: usize = 8;

/// BasesScorer reads four of the eight contigs at 70 completeness and any smaller group under 68.
fn bar() -> Rung {
    Rung {
        floor: 1,
        completeness: 68.0,
        contamination: 100.0,
    }
}

fn layout() -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    (
        Array2::zeros((CONTIGS, 2)),
        Array2::zeros((CONTIGS, 2)),
        vec![PIECE; CONTIGS],
    )
}

#[test]
fn a_bin_broken_into_pieces_that_all_miss_comes_back_whole() {
    let (coverage, tnf, lengths) = layout();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let dissolved = vec![(0usize, vec![0, 1, 2, 3])];

    let held = restore(
        &features,
        &BasesScorer::new(lengths.clone()),
        &dissolved,
        vec![vec![0, 1], vec![2, 3]],
        bar(),
    );

    assert_eq!(held.bins, 1);
    assert!(held.promoted.is_empty());
    let mut released = held.released;
    released.sort_unstable();
    assert_eq!(released, vec![0, 1, 2, 3]);
}

#[test]
fn a_bin_with_one_piece_over_the_bar_is_left_alone() {
    let (coverage, tnf, lengths) = layout();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let dissolved = vec![(0usize, vec![0, 1, 2, 3, 4, 5])];
    let promoted = vec![vec![0, 1, 2, 3], vec![4, 5]];

    let held = restore(
        &features,
        &BasesScorer::new(lengths.clone()),
        &dissolved,
        promoted.clone(),
        bar(),
    );

    assert_eq!(held.bins, 0);
    assert_eq!(held.promoted, promoted);
    assert!(held.released.is_empty());
}

/// Unpicking a piece that drew from two bins would orphan the other bin's contigs, so the
/// restore has to decline even though every test on the first bin passes.
#[test]
fn a_piece_drawing_from_two_bins_blocks_the_restore() {
    let (coverage, tnf, lengths) = layout();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let dissolved = vec![(0usize, vec![0, 1, 2, 3]), (1usize, vec![4, 5, 6, 7])];
    let promoted = vec![vec![0, 4]];

    let held = restore(
        &features,
        &BasesScorer::new(lengths.clone()),
        &dissolved,
        promoted.clone(),
        bar(),
    );

    assert_eq!(held.bins, 0);
    assert_eq!(held.promoted, promoted);
}

/// A dropped piece hands back everything it held, including the contigs it took from the
/// unbinned half of the pool, or they belong to no bin and no pool.
#[test]
fn a_dropped_piece_releases_the_contigs_it_took_from_the_pool() {
    let (coverage, tnf, lengths) = layout();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let dissolved = vec![(0usize, vec![0, 1, 2, 3])];

    let held = restore(
        &features,
        &BasesScorer::new(lengths.clone()),
        &dissolved,
        vec![vec![0, 1, 6]],
        bar(),
    );

    assert_eq!(held.bins, 1);
    let mut released = held.released;
    released.sort_unstable();
    assert_eq!(released, vec![0, 1, 6]);
}
