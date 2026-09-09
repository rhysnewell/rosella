use std::collections::BTreeMap;

use crate::embedding::features::ContigFeatures;

/// Two strains at one depth look alike in every feature, and so do two halves of one genome.
/// Only scale separates them, and the run's own closed genomes say what genome-sized is here.
pub fn floor(
    features: &ContigFeatures,
    bins: &BTreeMap<usize, Vec<usize>>,
    unbinned: &[usize],
    min_bin_size: usize,
) -> Option<usize> {
    let mut alone = bins
        .values()
        .filter(|&contigs| contigs.len() == 1 ).map(|contigs| features.length(contigs[0]))
        .chain(unbinned.iter().map(|contig| features.length(*contig)))
        .filter(|length| *length >= min_bin_size)
        .collect::<Vec<_>>();
    if alone.is_empty() {
        return None;
    }
    alone.sort_unstable();
    Some(alone[alone.len() / 2] / 2)
}
