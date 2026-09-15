//! Reverting a run of dissolved bins the pool made nothing better out of.

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::quality::Worth;
use rosella::refine::restore::restore;
use rosella::refine::rung::Rung;

#[path = "../support/scorer.rs"]
mod scorer;

use scorer::GenomeScorer;

const PIECE: usize = 100_000;
const CONTIGS: usize = 8;

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

fn bar() -> Rung {
    Rung {
        floor: 1,
        completeness: 90.0,
        contamination: 5.0,
    }
}

fn layout() -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    (
        Array2::zeros((CONTIGS, 2)),
        Array2::zeros((CONTIGS, 2)),
        vec![PIECE; CONTIGS],
    )
}

macro_rules! held {
    ($dissolved:expr, $promoted:expr) => {{
        let (coverage, tnf, lengths) = layout();
        let features = ContigFeatures::new(&coverage, &tnf, &lengths);
        restore(&features, &truth(), worth(), bar(), &$dissolved, $promoted)
    }};
}

#[test]
fn a_whole_genome_cut_in_half_goes_back_together() {
    let held = held!(vec![(0usize, vec![0, 1, 2, 3])], vec![vec![0, 1], vec![2, 3]]);

    assert_eq!(held.bins, 1);
    assert!(held.promoted.is_empty());
    let mut released = held.released;
    released.sort_unstable();
    assert_eq!(released, vec![0, 1, 2, 3]);
}

/// Splitting a fused bin is what the pool is for, and the two halves are genomes where the bin
/// was neither, so nothing is put back.
#[test]
fn a_fused_bin_split_into_two_genomes_is_left_alone() {
    let promoted = vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7]];
    let held = held!(vec![(0usize, vec![0, 1, 2, 3, 4, 5, 6, 7])], promoted.clone());

    assert_eq!(held.bins, 0);
    assert_eq!(held.promoted, promoted);
}

/// One piece drawing from two bins ties them into a single run, and reverting the run is worth
/// two genomes against the one bad bin the pool built.
#[test]
fn a_piece_spanning_two_bins_reverts_both_of_them() {
    let dissolved = vec![(0usize, vec![0, 1, 2, 3]), (1usize, vec![4, 5, 6, 7])];
    let held = held!(dissolved, vec![vec![0, 4]]);

    assert_eq!(held.bins, 2);
    assert!(held.promoted.is_empty());
    let mut released = held.released;
    released.sort_unstable();
    assert_eq!(released, vec![0, 4]);
}

/// The same run of two bins, but here the pool cut the genomes out of them. Unravelling would
/// destroy two genomes to rebuild two fused bins, so the run stands.
#[test]
fn a_run_the_pool_sorted_into_genomes_stands() {
    let dissolved = vec![(0usize, vec![0, 1, 4, 5]), (1usize, vec![2, 3, 6, 7])];
    let promoted = vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7]];
    let held = held!(dissolved, promoted.clone());

    assert_eq!(held.bins, 0);
    assert_eq!(held.promoted, promoted);
}

/// A dropped piece hands back everything it held, including the contigs it took from the
/// unbinned half of the pool, or they belong to no bin and no pool.
#[test]
fn a_dropped_piece_releases_the_contigs_it_took_from_the_pool() {
    let held = held!(vec![(0usize, vec![0, 1, 2, 3])], vec![vec![0, 1, 6]]);

    assert_eq!(held.bins, 1);
    let mut released = held.released;
    released.sort_unstable();
    assert_eq!(released, vec![0, 1, 6]);
}
