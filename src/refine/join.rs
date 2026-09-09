use std::collections::BTreeMap;

use crate::embedding::features::ContigFeatures;
use crate::quality::ContigQuality;

fn reciprocated(best: &[Option<(f64, usize)>]) -> Vec<(f64, usize, usize)> {
    let mut pairs = Vec::new();
    for (left, nearest) in best.iter().enumerate() {
        let Some((distance, right)) = *nearest else {
            continue;
        };
        if left < right && best[right].is_some_and(|(_, back)| back == left) {
            pairs.push((distance, left, right));
        }
    }
    pairs
}

/// A pair joins, then the pair it made can take a third piece, but the chain is short and every
/// pass costs a full sweep of the boosters.
const MAX_PASSES: usize = 4;

#[derive(Debug, Default, Clone, Copy)]
pub struct JoinLedger {
    pub bins: usize,
    pub short: usize,
    pub scored: usize,
    pub reciprocated: usize,
    pub joined: usize,
    pub passes: usize,
}

impl std::fmt::Display for JoinLedger {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "{} bins over {} passes, {} short of whole; scored {} unions, {} pairs each \
             other's best, joined {}",
            self.bins, self.passes, self.short, self.scored, self.reciprocated, self.joined
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

struct Piece {
    contigs: Vec<usize>,
    completeness: f64,
    families: std::collections::HashSet<u32>,
}

/// Two halves of one genome hold different gene families, so their union is more complete than
/// either and no more contaminated. Neither composition nor coverage separates such a pair from
/// the far larger number of pairs that merely sit close together.
fn pass(
    features: &ContigFeatures,
    quality: &ContigQuality,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    settings: JoinSettings,
    ledger: &mut JoinLedger,
) -> usize {
    // A bin the model already calls whole has nothing to gain and everything to lose, so it
    // neither takes a partner nor is offered as one.
    let mut ids = Vec::new();
    let mut pieces = Vec::new();
    for id in bins.keys().copied() {
        let contigs = bins[&id].clone();
        let completeness = quality.score(&contigs).completeness;
        if completeness >= settings.completeness {
            continue;
        }
        ids.push(id);
        pieces.push(Piece {
            completeness,
            families: quality.families(&contigs),
            contigs,
        });
    }

    ledger.short = ledger.short.max(pieces.len());
    let mut best = vec![None; pieces.len()];
    for left in 0..pieces.len() {
        for right in (left + 1)..pieces.len() {
            let novel = pieces[right]
                .families
                .difference(&pieces[left].families)
                .count();
            let lacking = pieces[left]
                .families
                .difference(&pieces[right].families)
                .count();
            if novel == 0 || lacking == 0 {
                continue;
            }
            let joined = union(&pieces[left].contigs, &pieces[right].contigs);
            if features.bin_size(&joined) > settings.max_bin_size {
                continue;
            }
            let held = quality.score(&joined);
            ledger.scored += 1;
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
    }

    let mut pairs = reciprocated(&best);
    ledger.reciprocated += pairs.len();
    pairs.sort_by(|one, other| other.0.total_cmp(&one.0));

    let mut taken = vec![false; pieces.len()];
    let mut joined = 0;
    for (_, left, right) in pairs {
        if taken[left] || taken[right] {
            continue;
        }
        taken[right] = true;
        let contigs = union(&pieces[left].contigs, &pieces[right].contigs);
        pieces[left].contigs = contigs.clone();
        bins.insert(ids[left], contigs);
        bins.remove(&ids[right]);
        joined += 1;
    }
    joined
}

pub fn join(
    features: &ContigFeatures,
    quality: &ContigQuality,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    settings: JoinSettings,
) -> JoinLedger {
    let mut ledger = JoinLedger::default();
    for _ in 0..MAX_PASSES {
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
