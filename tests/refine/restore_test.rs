//! Putting a dissolved bin back when the pool made nothing better out of it.

use rosella::quality::Worth;
use rosella::refine::restore::restore;

#[path = "../support/scorer.rs"]
mod scorer;

use scorer::GenomeScorer;

/// Contigs 0 to 3 are one genome, 4 to 7 another.
fn truth() -> GenomeScorer {
    GenomeScorer::new(vec![0, 0, 0, 0, 1, 1, 1, 1])
}

fn worth() -> Worth {
    Worth {
        contamination: 2.0,
        allowance: 0.0,
    }
}

#[test]
fn a_whole_genome_cut_in_half_goes_back_together() {
    let dissolved = vec![(0usize, vec![0, 1, 2, 3])];

    let held = restore(
        &truth(),
        worth(),
        &dissolved,
        vec![vec![0, 1], vec![2, 3]],
    );

    assert_eq!(held.bins, 1);
    assert!(held.promoted.is_empty());
    let mut released = held.released;
    released.sort_unstable();
    assert_eq!(released, vec![0, 1, 2, 3]);
}

/// Splitting a fused bin is what the pool is for, and both halves outscore the bin they came
/// from, so nothing is put back.
#[test]
fn a_fused_bin_split_into_two_genomes_is_left_alone() {
    let dissolved = vec![(0usize, vec![0, 1, 2, 3, 4, 5, 6, 7])];
    let promoted = vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7]];

    let held = restore(&truth(), worth(), &dissolved, promoted.clone());

    assert_eq!(held.bins, 0);
    assert_eq!(held.promoted, promoted);
    assert!(held.released.is_empty());
}

/// Unpicking a piece that drew from two bins would orphan the other bin's contigs, so the
/// restore declines even though the first bin was worth more whole.
#[test]
fn a_piece_drawing_from_two_bins_blocks_the_restore() {
    let dissolved = vec![(0usize, vec![0, 1, 2, 3]), (1usize, vec![4, 5, 6, 7])];
    let promoted = vec![vec![0, 4]];

    let held = restore(&truth(), worth(), &dissolved, promoted.clone());

    assert_eq!(held.bins, 0);
    assert_eq!(held.promoted, promoted);
}

/// A dropped piece hands back everything it held, including the contigs it took from the
/// unbinned half of the pool, or they belong to no bin and no pool.
#[test]
fn a_dropped_piece_releases_the_contigs_it_took_from_the_pool() {
    let dissolved = vec![(0usize, vec![0, 1, 2, 3])];

    let held = restore(&truth(), worth(), &dissolved, vec![vec![0, 1, 6]]);

    assert_eq!(held.bins, 1);
    let mut released = held.released;
    released.sort_unstable();
    assert_eq!(released, vec![0, 1, 6]);
}
