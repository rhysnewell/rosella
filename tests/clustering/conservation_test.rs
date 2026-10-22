//! Every stage partitions the contigs exactly and nothing checked it, so this is the check.

use std::collections::HashSet;

use rosella::clustering::clusterer::{conserved, placed_once};

fn four() -> HashSet<usize> {
    (0..4).collect()
}

#[test]
fn a_contig_placed_twice_or_from_nowhere_is_caught() {
    assert!(placed_once([0, 1, 1], &four()).is_err());
    assert!(placed_once([0, 9], &four()).is_err());
    assert_eq!(placed_once([0, 1], &four()).unwrap().len(), 2);
}

#[test]
fn a_contig_that_never_came_out_of_binning_is_caught() {
    assert!(conserved([0, 1, 2], &four()).is_err());
    assert!(conserved([0, 1, 2, 3], &four()).is_ok());
}
