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

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum Gate {
    #[default]
    Bars,
    Excess,
    Sized,
}

#[derive(Default)]
struct Tally {
    short: u32,
    adrift: u32,
    short_bp: usize,
    adrift_bp: usize,
    bp: usize,
    anchored: bool,
}

// A bin that already clears the bars gains nothing by losing a contig, and its own short contigs
// are the ones its completeness rests on.
pub fn adrift(
    bins: &BTreeMap<usize, Vec<usize>>,
    finished: &HashSet<usize>,
    lengths: &[usize],
    knn: &KnnGraph,
    gate: Gate,
    contamination: f64,
) -> Vec<usize> {
    let owner = owners(bins.iter().map(|(bin, members)| (*bin, members)));
    let mut adrift = owner
        .par_iter()
        .filter(|(contig, bin)| {
            lengths[**contig] < LEAST_BASES && (gate != Gate::Bars || !finished.contains(*bin)) && {
                let mut weights = neighbour_weight(**contig, &owner, knn, lengths);
                share(&mut weights, **bin).is_some_and(|share| share < ADRIFT)
            }
        })
        .map(|(contig, _)| *contig)
        .collect::<Vec<_>>();
    adrift.sort_unstable();
    if gate == Gate::Bars {
        return adrift;
    }
    let shed = shedding(bins, finished, lengths, &adrift, gate, contamination);
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
    gate: Gate,
    contamination: f64,
) -> HashSet<usize> {
    let tallies = bins
        .iter()
        .map(|(bin, members)| {
            let mut tally = Tally::default();
            for contig in members {
                let length = lengths[*contig];
                tally.bp += length;
                if length >= LEAST_BASES {
                    tally.anchored = true;
                    continue;
                }
                tally.short += 1;
                tally.short_bp += length;
                if adrift.binary_search(contig).is_ok() {
                    tally.adrift += 1;
                    tally.adrift_bp += length;
                }
            }
            (*bin, tally)
        })
        .collect::<BTreeMap<_, _>>();
    let clean = tallies
        .iter()
        .filter(|(bin, _)| finished.contains(*bin))
        .map(|(_, tally)| tally);
    let (adrift_n, short_n, adrift_bp, short_bp) = clean.fold((0, 0, 0, 0), |sum, tally| {
        (
            sum.0 + tally.adrift,
            sum.1 + tally.short,
            sum.2 + tally.adrift_bp,
            sum.3 + tally.short_bp,
        )
    });
    if short_n == 0 {
        return HashSet::new();
    }
    let rate = f64::from(adrift_n) / f64::from(short_n);
    let rate_bp = adrift_bp as f64 / short_bp as f64;
    let open = tallies
        .iter()
        .filter(|(bin, tally)| tally.anchored && !finished.contains(*bin))
        .collect::<Vec<_>>();
    let tested = open.len().max(1) as f64;
    let shed = open
        .into_iter()
        .filter(|(_, tally)| {
            tally.adrift > 0
                && upper_tail(tally.adrift, tally.short, rate) < FAMILY_ERROR / tested
                && (gate != Gate::Sized
                    || tally.adrift_bp as f64 - rate_bp * tally.short_bp as f64
                        > contamination / 100.0 * tally.bp as f64)
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
