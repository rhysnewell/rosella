use std::collections::{HashMap, HashSet};

use crate::markers::ContigMarkers;
use crate::quality::{Bars, Scorer};

/// A bin that holds two whole copies of a single copy marker is holding sequence from two
/// genomes, and the copy that brings nothing else is the one that can leave. A bin already over
/// both bars has no second genome the markers can see, so thinning it only costs it sequence.
pub fn shed(
    bins: &mut HashMap<usize, HashSet<usize>>,
    unbinned: &mut HashSet<usize>,
    markers: &ContigMarkers,
    bars: Bars,
) -> usize {
    let mut labels = bins.keys().copied().collect::<Vec<_>>();
    labels.sort_unstable();
    let mut dropped = 0;
    for label in labels {
        let Some(members) = bins.get_mut(&label) else {
            continue;
        };
        let mut contigs = members.iter().copied().collect::<Vec<_>>();
        contigs.sort_unstable();
        if markers.score(&contigs).clears(bars) {
            continue;
        }
        for contig in markers.redundant(&contigs) {
            members.remove(&contig);
            unbinned.insert(contig);
            dropped += 1;
        }
    }
    bins.retain(|_, members| !members.is_empty());
    dropped
}
