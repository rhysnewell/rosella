use std::collections::BTreeMap;

use crate::markers::{ContigMarkers, Strict};
use crate::quality::{Bars, Scorer};

// Two whole copies of a single copy marker mean two genomes, and the copy that brings nothing
// else is the one that does not belong.
pub fn shed(
    bins: &mut BTreeMap<usize, Vec<usize>>,
    unbinned: &mut Vec<usize>,
    markers: &ContigMarkers,
    bars: Bars,
    strict: Option<&Strict>,
) -> usize {
    let mut dropped = 0;
    for members in bins.values_mut() {
        members.sort_unstable();
        let plain = markers.score(members);
        let judged = strict.and_then(|strict| strict.judged(members, plain));
        if judged.unwrap_or(plain).clears(bars) {
            continue;
        }
        let mut redundant = markers.redundant_on(members, judged.is_some());
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
