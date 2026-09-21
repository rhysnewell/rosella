use std::collections::BTreeMap;

use crate::embedding::features::ContigFeatures;
use crate::markers::ContigMarkers;
use crate::quality::{Bars, Scorer};
use crate::refine::bisect;

pub struct Split<'a> {
    pub features: &'a ContigFeatures<'a>,
    pub min_bin_size: usize,
    pub seed: u64,
}

/// A bin that holds two whole copies of a single copy marker is holding sequence from two
/// genomes, and the copy that brings nothing else is the one that can leave. A bin already over
/// both bars has no second genome the markers can see, so thinning it only costs it sequence.
pub fn shed(
    bins: &mut BTreeMap<usize, Vec<usize>>,
    unbinned: &mut Vec<usize>,
    markers: &ContigMarkers,
    bars: Bars,
    spacings: f64,
    held: &dyn Fn(&[usize]) -> bool,
    split: Option<Split<'_>>,
) -> usize {
    let skip = |members: &[usize]| markers.score(members).clears(bars) || held(members);
    let fused = match split {
        Some(_) => bins.values().filter(|members| !skip(members)).count(),
        None => 0,
    };
    let mut next = bins.keys().copied().max().map_or(0, |label| label + 1);
    let mut grown = Vec::new();
    let mut dropped = 0;
    for members in bins.values_mut() {
        members.sort_unstable();
        if skip(members) {
            continue;
        }
        if let Some(split) = &split
            && let Some([kept, moved]) = partition(split, markers, members, fused, spacings)
        {
            *members = kept;
            grown.push(moved);
            continue;
        }
        let mut redundant = markers.redundant(members, spacings);
        if redundant.is_empty() {
            continue;
        }
        redundant.sort_unstable();
        let kept = members
            .iter()
            .copied()
            .filter(|contig| redundant.binary_search(contig).is_err())
            .collect::<Vec<_>>();
        *members = kept;
        dropped += redundant.len();
        unbinned.append(&mut redundant);
    }
    for piece in grown {
        bins.insert(next, piece);
        next += 1;
    }
    bins.retain(|_, members| !members.is_empty());
    dropped
}

/// The victim and the carrier that made it look redundant are the markers' own nomination of
/// the two genomes, so the split grows from them rather than from the geometry's extremes.
fn partition(
    split: &Split<'_>,
    markers: &ContigMarkers,
    members: &[usize],
    fused: usize,
    spacings: f64,
) -> Option<[Vec<usize>; 2]> {
    let first = markers.redundant_traced(members, spacings).into_iter().next()?;
    let twin = first.twin?;
    let victim = members.iter().position(|contig| *contig == first.contig)?;
    let carrier = members.iter().position(|contig| *contig == twin)?;
    bisect::from_seeds(
        split.features,
        members,
        (victim, carrier),
        split.min_bin_size,
        fused,
        split.seed,
    )
}
