use std::collections::{BTreeMap, HashMap, HashSet};

use rayon::prelude::*;

use crate::embedding::knn::KnnGraph;
use crate::markers::replicon::LEAST_BASES;
use crate::quality::{FAMILY_ERROR, upper_tail};
use crate::refine::neighbourhood::{neighbour_weight, share};
use crate::refine::owners::owners;

// The partition labels every contig, so a short one whose neighbours mostly sit in other bins
// was placed on weak evidence.
const ADRIFT: f64 = 0.5;

pub struct Reading {
    pub contig: usize,
    pub bin: usize,
    pub weighted: f64,
    pub counted: f64,
    pub top: f64,
}

#[derive(Default)]
struct Tally {
    short: u32,
    adrift: u32,
    anchored: bool,
}

// A bin that already clears the bars gains nothing by losing a contig, and its own short contigs
// are the ones its completeness rests on.
pub fn adrift(
    bins: &BTreeMap<usize, Vec<usize>>,
    finished: &HashSet<usize>,
    lengths: &[usize],
    knn: &KnnGraph,
) -> Vec<usize> {
    let owner = owners(bins.iter().map(|(bin, members)| (*bin, members)));
    let mut adrift = owner
        .par_iter()
        .filter(|(contig, bin)| {
            lengths[**contig] < LEAST_BASES && {
                let mut weights = neighbour_weight(**contig, &owner, knn, lengths);
                share(&mut weights, **bin).is_some_and(|share| share < ADRIFT)
            }
        })
        .map(|(contig, _)| *contig)
        .collect::<Vec<_>>();
    adrift.sort_unstable();
    let shed = shedding(bins, finished, lengths, &adrift);
    adrift.retain(|contig| shed.contains(&owner[contig]));
    adrift
}

// Every bin holds a few odd pieces of its own that sit nearer other bins, at the rate the finished
// bins show. Only a bin holding more than that rate allows sheds them.
fn shedding(
    bins: &BTreeMap<usize, Vec<usize>>,
    finished: &HashSet<usize>,
    lengths: &[usize],
    adrift: &[usize],
) -> HashSet<usize> {
    let tallies = bins
        .iter()
        .map(|(bin, members)| {
            let mut tally = Tally::default();
            for contig in members {
                if lengths[*contig] >= LEAST_BASES {
                    tally.anchored = true;
                    continue;
                }
                tally.short += 1;
                tally.adrift += u32::from(adrift.binary_search(contig).is_ok());
            }
            (*bin, tally)
        })
        .collect::<BTreeMap<_, _>>();
    let (adrift_n, short_n) = tallies
        .iter()
        .filter(|(bin, _)| finished.contains(*bin))
        .fold((0, 0), |sum, (_, tally)| {
            (sum.0 + tally.adrift, sum.1 + tally.short)
        });
    if short_n == 0 {
        return HashSet::new();
    }
    let rate = f64::from(adrift_n) / f64::from(short_n);
    let open = tallies
        .iter()
        .filter(|(bin, tally)| tally.anchored && !finished.contains(*bin))
        .collect::<Vec<_>>();
    let tested = open.len().max(1) as f64;
    open.into_iter()
        .filter(|(_, tally)| {
            tally.adrift > 0 && upper_tail(tally.adrift, tally.short, rate) < FAMILY_ERROR / tested
        })
        .map(|(bin, _)| *bin)
        .collect::<HashSet<_>>()
}

pub fn readings(
    bins: &BTreeMap<usize, Vec<usize>>,
    lengths: &[usize],
    knn: &KnnGraph,
) -> Vec<Reading> {
    let owner = owners(bins.iter().map(|(bin, members)| (*bin, members)));
    let mut out = owner
        .par_iter()
        .filter(|(contig, _)| lengths[**contig] < LEAST_BASES && **contig < knn.n_points())
        .filter_map(|(contig, bin)| {
            let mut weights = neighbour_weight(*contig, &owner, knn, lengths);
            let weighted = share(&mut weights, *bin)?;
            let total = weights.iter().map(|(_, weight)| *weight).sum::<f64>();
            let top = weights
                .iter()
                .filter(|(other, _)| other != bin)
                .map(|(_, weight)| *weight / total)
                .fold(0.0, f64::max);
            let held = knn
                .indices
                .row(*contig)
                .iter()
                .filter_map(|neighbour| owner.get(&(*neighbour as usize)))
                .collect::<Vec<_>>();
            let counted = held.iter().filter(|label| **label == bin).count() as f64
                / held.len().max(1) as f64;
            Some(Reading {
                contig: *contig,
                bin: *bin,
                weighted,
                counted,
                top,
            })
        })
        .collect::<Vec<_>>();
    out.sort_unstable_by_key(|reading| reading.contig);
    out
}

// A bin's own contigs can fill at most its size less one of a contig's neighbours, so a small
// bin whose contigs do not reach half of that was grouped on weak evidence.
pub fn scattered(members: &[usize], owner: &HashMap<usize, usize>, knn: &KnnGraph) -> bool {
    let Some(bin) = members.first().and_then(|contig| owner.get(contig)) else {
        return false;
    };
    let reach = members
        .iter()
        .filter(|contig| **contig < knn.n_points())
        .filter_map(|contig| {
            let held = knn
                .indices
                .row(*contig)
                .iter()
                .filter(|neighbour| **neighbour as usize != *contig)
                .filter_map(|neighbour| owner.get(&(*neighbour as usize)))
                .collect::<Vec<_>>();
            let room = held.len().min(members.len() - 1);
            (room > 0)
                .then(|| held.iter().filter(|label| **label == bin).count() as f64 / room as f64)
        })
        .collect::<Vec<_>>();
    !reach.is_empty() && reach.iter().sum::<f64>() / (reach.len() as f64) < ADRIFT
}
