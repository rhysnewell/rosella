//! Whether a fragment can rejoin the bin that already holds every marker family it carries.

use std::collections::BTreeMap;

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::join::{JoinSettings, Novelty, join};

#[path = "../support/scorer.rs"]
mod scorer;

use scorer::GenomeScorer;

const PIECE: usize = 100_000;

fn settings(novelty: Novelty) -> JoinSettings {
    JoinSettings {
        completeness: 90.0,
        contamination: 5.0,
        max_bin_size: 15_000_000,
        novelty,
    }
}

/// Eight contigs of one genome cut into a marker-rich bin and a marker-poor fragment whose
/// families the rich bin already carries. Cutting a genome in two is exactly how the pool
/// makes this pair, and the mutual novelty test is what refuses it.
fn cut_genome() -> (Vec<usize>, Vec<u32>, BTreeMap<usize, Vec<usize>>) {
    let genomes = vec![0; 8];
    let families = vec![0, 1, 2, 3, 4, 0, 1, 2];
    let bins = BTreeMap::from([(0usize, vec![0, 1, 2, 3, 4]), (1usize, vec![5, 6, 7])]);
    (genomes, families, bins)
}

#[test]
fn a_subset_fragment_rejoins_its_parent_only_once_novelty_is_off() {
    let coverage = Array2::zeros((8, 2));
    let tnf = Array2::zeros((8, 2));
    let lengths = vec![PIECE; 8];
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);

    let (genomes, families, strict) = cut_genome();
    let quality = GenomeScorer::new(genomes).with_families(families);

    let mut bins = strict.clone();
    join(&features, &quality, &mut bins, settings(Novelty::Strict));
    assert_eq!(bins, strict, "the strict test never scores the union");

    let mut bins = strict;
    let ledger = join(&features, &quality, &mut bins, settings(Novelty::Gain));
    assert_eq!(ledger.joined, 1);
    assert_eq!(bins.len(), 1);
    assert_eq!(bins[&0], (0..8).collect::<Vec<_>>());
}
