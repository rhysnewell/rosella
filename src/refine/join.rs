use std::collections::{BTreeMap, HashSet};

use rayon::prelude::*;

use crate::embedding::features::ContigFeatures;
use crate::quality::Scorer;

fn reciprocated(best: &[Option<(f64, usize)>]) -> Vec<(usize, usize)> {
    best.iter()
        .enumerate()
        .filter_map(|(left, nearest)| {
            let (_, right) = (*nearest)?;
            (left < right && best[right].is_some_and(|(_, back)| back == left))
                .then_some((left, right))
        })
        .collect()
}

#[derive(Debug, Default, Clone, Copy)]
pub struct JoinLedger {
    pub bins: usize,
    pub short: usize,
    pub scored: usize,
    pub joined: usize,
    pub passes: usize,
}

impl std::fmt::Display for JoinLedger {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "{} bins over {} passes, {} short of whole; scored {} unions, joined {} pairs that \
             were each other's best",
            self.bins, self.passes, self.short, self.scored, self.joined
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct JoinSettings {
    pub completeness: f64,
    pub contamination: f64,
    pub max_bin_size: usize,
}

fn union(left: &[usize], right: &[usize]) -> Vec<usize> {
    let mut joined = Vec::with_capacity(left.len() + right.len());
    joined.extend_from_slice(left);
    joined.extend_from_slice(right);
    joined.sort_unstable();
    joined
}

struct Piece<'a> {
    contigs: &'a [usize],
    bases: usize,
    completeness: f64,
    short: bool,
    families: HashSet<u32>,
}

// Two halves of one genome hold different markers, so their union is more complete and no more
// contaminated. Neither composition nor coverage tells such a pair from merely close ones.
fn pass(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    settings: JoinSettings,
    ledger: &mut JoinLedger,
) -> usize {
    // Marker completeness cannot see bases a bin is missing when that sequence carries no
    // marker, so a bin the model calls whole still enters the list as a receiver.
    let ids = bins.keys().copied().collect::<Vec<_>>();
    let pieces = bins
        .values()
        .map(|contigs| {
            let completeness = quality.score(contigs).completeness;
            Piece {
                completeness,
                short: completeness < settings.completeness,
                bases: features.bin_size(contigs),
                families: quality.features(contigs),
                contigs,
            }
        })
        .collect::<Vec<_>>();

    ledger.short = ledger
        .short
        .max(pieces.iter().filter(|piece| piece.short).count());
    let unions = (0..pieces.len())
        .into_par_iter()
        .flat_map_iter(|left| {
            let pieces = &pieces;
            ((left + 1)..pieces.len()).filter_map(move |right| {
                let (one, other) = (&pieces[left], &pieces[right]);
                if !one.short && !other.short {
                    return None;
                }
                if one.bases + other.bases > settings.max_bin_size {
                    return None;
                }
                // Refuses half the unions. On all 20 measured sets the bars below reject every
                // one of them anyway, so this only saves the scoring.
                if other.families.is_subset(&one.families)
                    || one.families.is_subset(&other.families)
                {
                    return None;
                }
                let held = quality.score(&union(one.contigs, other.contigs));
                Some((left, right, held))
            })
        })
        .collect::<Vec<_>>();
    ledger.scored += unions.len();

    let mut best = vec![None; pieces.len()];
    for (left, right, held) in unions {
        if held.contamination > settings.contamination {
            continue;
        }
        let gain = held.completeness - pieces[left].completeness.max(pieces[right].completeness);
        if gain <= 0.0 {
            continue;
        }
        if best[left].is_none_or(|(most, _)| gain > most) {
            best[left] = Some((gain, right));
        }
        if best[right].is_none_or(|(most, _)| gain > most) {
            best[right] = Some((gain, left));
        }
    }

    // Every piece has one best, so no piece sits in two reciprocated pairs.
    let joins = reciprocated(&best)
        .into_iter()
        .map(|(left, right)| {
            let contigs = union(pieces[left].contigs, pieces[right].contigs);
            (ids[left], ids[right], contigs)
        })
        .collect::<Vec<_>>();
    let joined = joins.len();
    for (left, right, contigs) in joins {
        bins.insert(left, contigs);
        bins.remove(&right);
    }
    joined
}

pub fn join(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    settings: JoinSettings,
) -> JoinLedger {
    let mut ledger = JoinLedger::default();
    for _ in 0..crate::tuning::JOIN_PASSES {
        ledger.passes += 1;
        let joined = pass(features, quality, bins, settings, &mut ledger);
        ledger.joined += joined;
        if joined == 0 {
            break;
        }
    }
    ledger.bins = bins.len();
    ledger
}
