use std::collections::BTreeMap;

use crate::markers::ContigMarkers;
use crate::quality::{Bars, Scorer};

/// A bin that holds two whole copies of a single copy marker is holding sequence from two
/// genomes, and the copy that brings nothing else is the one that can leave. A bin already over
/// both bars has no second genome the markers can see, so thinning it only costs it sequence.
pub fn shed(
    bins: &mut BTreeMap<usize, Vec<usize>>,
    unbinned: &mut Vec<usize>,
    markers: &ContigMarkers,
    bars: Bars,
) -> usize {
    let mut dropped = 0;
    for members in bins.values_mut() {
        members.sort_unstable();
        if markers.score(members).clears(bars) {
            continue;
        }
        let mut redundant = markers.redundant(members);
        if redundant.is_empty() {
            continue;
        }
        redundant.sort_unstable();
        members.retain(|contig| redundant.binary_search(contig).is_err());
        dropped += redundant.len();
        unbinned.append(&mut redundant);
    }
    bins.retain(|_, members| !members.is_empty());
    dropped
}
