use std::collections::BTreeMap;

use log::debug;

use crate::embedding::features::ContigFeatures;

/// Two strains at one depth look alike in every feature, and so do two halves of one genome.
/// Only scale separates them, and the run's own closed genomes say what genome-sized is here.
pub fn floor(
    features: &ContigFeatures,
    bins: &BTreeMap<usize, Vec<usize>>,
    unbinned: &[usize],
    min_bin_size: usize,
) -> Option<usize> {
    let alone = bins
        .values()
        .filter(|&contigs| contigs.len() == 1)
        .map(|contigs| features.length(contigs[0]))
        .chain(unbinned.iter().map(|contig| features.length(*contig)))
        .filter(|length| *length >= min_bin_size)
        .collect::<Vec<_>>();
    debug!("Genome scale read off {} contigs binned alone", alone.len());
    median(alone).map(|middle| middle / 2)
}


fn median(mut values: Vec<usize>) -> Option<usize> {
    if values.is_empty() {
        return None;
    }
    values.sort_unstable();
    Some(values[values.len() / 2])
}
