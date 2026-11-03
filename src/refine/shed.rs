use std::collections::BTreeMap;

use crate::markers::{ContigMarkers, Strict};
use crate::quality::{Bars, Scorer};

// Two whole copies of a single copy marker mean two genomes, and the copy that brings nothing
// else can leave. A bin over both bars has no second genome to see, so thinning it only costs.
pub fn shed(
    bins: &mut BTreeMap<usize, Vec<usize>>,
    unbinned: &mut Vec<usize>,
    markers: &ContigMarkers,
    bars: Bars,
    checkm: bool,
) -> usize {
    let strict = Strict {
        inner: markers,
        markers,
    };
    let scorer: &dyn Scorer = match checkm {
        true => &strict,
        false => markers,
    };
    let mut dropped = 0;
    for members in bins.values_mut() {
        members.sort_unstable();
        if scorer.score(members).clears(bars) {
            continue;
        }
        let mut redundant = markers.redundant_on(members, checkm);
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
