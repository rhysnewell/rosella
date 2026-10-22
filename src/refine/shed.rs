use std::collections::BTreeMap;

use crate::markers::ContigMarkers;
use crate::quality::{Bars, Scorer};

// Two whole copies of a single copy marker mean two genomes, and the copy that brings nothing
// else can leave. A bin over both bars has no second genome to see, so thinning it only costs.
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
        let kept = members
            .iter()
            .copied()
            .filter(|contig| redundant.binary_search(contig).is_err())
            .collect::<Vec<_>>();
        *members = kept;
        dropped += redundant.len();
        unbinned.append(&mut redundant);
    }
    bins.retain(|_, members| !members.is_empty());
    dropped
}
