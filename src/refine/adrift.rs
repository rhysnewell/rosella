use std::collections::{BTreeMap, HashSet};

use log::info;
use rayon::prelude::*;

use crate::embedding::knn::KnnGraph;
use crate::markers::replicon::LEAST_BASES;
use crate::quality::{FAMILY_ERROR, upper_tail};
use crate::refine::owners::owners;
use crate::refine::report_context::{neighbour_weight, share};

// The partition labels every contig, so a short one whose neighbours mostly sit in other bins
// was placed on weak evidence.
const ADRIFT: f64 = 0.5;

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
    let shed = open
        .into_iter()
        .filter(|(_, tally)| {
            tally.adrift > 0 && upper_tail(tally.adrift, tally.short, rate) < FAMILY_ERROR / tested
        })
        .map(|(bin, _)| *bin)
        .collect::<HashSet<_>>();
    info!(
        "Adrift: {:.3} of short contigs in finished bins, {} of {} open bins shed theirs",
        rate,
        shed.len(),
        tested
    );
    shed
}
