use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};

use anyhow::Result;
use log::warn;

use crate::clustering::clusterer::Partitioning;
use crate::refine::dissolve::{DissolveLedger, DissolveSettings, Pot, RoundParams, neighbours_for};
use crate::refine::rung::{RUNGS, Rung, Verdict};

struct Ranked {
    worth: f64,
    contigs: Vec<usize>,
    verdict: Verdict,
}

impl PartialEq for Ranked {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Ranked {}

impl PartialOrd for Ranked {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Ranked {
    fn cmp(&self, other: &Self) -> Ordering {
        self.worth
            .total_cmp(&other.worth)
            .then_with(|| other.contigs.cmp(&self.contigs))
    }
}

/// A proposal the winners already emptied is not the proposal that was scored, so what is left
/// of it is put back through the same bar rather than trusted on the rank it earned whole.
pub fn remaining(contigs: &[usize], claimed: &HashSet<usize>) -> Vec<usize> {
    contigs
        .iter()
        .copied()
        .filter(|contig| !claimed.contains(contig))
        .collect()
}

fn remaining_in(contigs: &[usize], pool: &HashSet<usize>) -> Vec<usize> {
    contigs
        .iter()
        .copied()
        .filter(|contig| pool.contains(contig))
        .collect()
}

fn sorted(contigs: HashSet<usize>) -> Vec<usize> {
    let mut contigs = contigs.into_iter().collect::<Vec<_>>();
    contigs.sort_unstable();
    contigs
}

fn propose(
    pool: &HashSet<usize>,
    settings: DissolveSettings,
    oracle: &[Vec<usize>],
    ledger: &mut DissolveLedger,
    partition: &impl Fn(&HashSet<usize>, RoundParams) -> Result<Vec<Partitioning>>,
) -> Vec<Vec<usize>> {
    let mut candidates = oracle
        .iter()
        .map(|group| remaining_in(group, pool))
        .filter(|group| group.len() >= 2)
        .collect::<Vec<_>>();
    for round in 0..settings.rounds.max(1) {
        let results = match partition(pool, neighbours_for(settings, round)) {
            Ok(results) => results,
            Err(error) => {
                warn!("Could not re-embed the pool: {error}");
                break;
            }
        };
        ledger.rounds += 1;
        for result in results {
            ledger.noise = result.outliers.len();
            candidates.extend(result.cluster_map.into_values().map(sorted));
        }
    }
    candidates
}

fn dedupe(candidates: &mut Vec<Vec<usize>>) {
    candidates.sort_unstable();
    candidates.dedup();
}

fn heap(pot: &Pot, candidates: Vec<Vec<usize>>) -> BinaryHeap<Ranked> {
    candidates
        .into_iter()
        .map(|contigs| Ranked {
            worth: pot.worth(&contigs),
            contigs,
            verdict: Verdict::Adopt,
        })
        .collect()
}

fn sweep(
    pot: &Pot,
    mut held: BinaryHeap<Ranked>,
    claimed: &mut HashSet<usize>,
    bar: Rung,
) -> (Vec<Vec<usize>>, BinaryHeap<Ranked>) {
    let mut taken = Vec::new();
    let mut refused = BinaryHeap::new();
    while let Some(entry) = held.pop() {
        let left = remaining(&entry.contigs, claimed);
        if left.len() < 2 {
            continue;
        }
        let verdict = match pot.judge(&left, bar) {
            Verdict::Adopt if !pot.improves(&left) => Verdict::Adopt,
            Verdict::Adopt => {
                claimed.extend(left.iter().copied());
                taken.push(left);
                continue;
            }
            other => other,
        };
        refused.push(Ranked {
            worth: entry.worth,
            contigs: left,
            verdict,
        });
    }
    (taken, refused)
}

fn tally(refused: &BinaryHeap<Ranked>, ledger: &mut DissolveLedger) {
    for entry in refused {
        match entry.verdict {
            Verdict::TooSmall => ledger.refused_small += 1,
            Verdict::Incomplete => ledger.refused_incomplete += 1,
            Verdict::Contaminated => ledger.refused_contaminated += 1,
            Verdict::Duplicated => ledger.refused_duplicated += 1,
            Verdict::Adopt => ledger.refused_worse += 1,
        }
    }
}

fn claim(
    pot: &Pot,
    candidates: Vec<Vec<usize>>,
    settings: DissolveSettings,
    top: usize,
    ledger: &mut DissolveLedger,
) -> Vec<Vec<usize>> {
    let scored = pot.scored();
    let mut promoted = Vec::new();
    let mut claimed = HashSet::new();
    let mut held = heap(pot, candidates);

    for at in ledger.rung..RUNGS {
        ledger.rung = at;
        let bar = settings.bars.at(top, at, scored);
        let (taken, refused) = sweep(pot, held, &mut claimed, bar);
        let empty = taken.is_empty();
        promoted.extend(taken);
        held = refused;
        if !empty {
            break;
        }
    }
    tally(&held, ledger);
    promoted
}

/// Every round searches the same pool, so the bar is asked which proposal to keep rather than
/// which came first, and each pass then re-embeds what the pass before it left.
pub fn ranked(
    pot: &Pot,
    pool: &mut HashSet<usize>,
    settings: DissolveSettings,
    oracle: &[Vec<usize>],
    top: usize,
    ledger: &mut DissolveLedger,
    partition: impl Fn(&HashSet<usize>, RoundParams) -> Result<Vec<Partitioning>>,
) -> Vec<Vec<usize>> {
    let mut promoted = Vec::new();
    let mut before: Option<f64> = None;
    for _ in 0..settings.passes.max(1) {
        if pool.len() < settings.min_contigs {
            break;
        }
        let mut candidates = propose(pool, settings, oracle, ledger, &partition);
        dedupe(&mut candidates);
        ledger.proposed += candidates.len();

        let taken = claim(pot, candidates, settings, top, ledger);
        if taken.is_empty() {
            break;
        }
        for contigs in &taken {
            for contig in contigs {
                pool.remove(contig);
            }
        }
        ledger.passes += 1;
        // How many passes a pool is worth differs per assembly, and the bins a pass finds are
        // worth less than the last one's long before it finds none, which no tier can see.
        let held = taken.iter().map(|contigs| pot.worth(contigs)).sum::<f64>() / taken.len() as f64;
        promoted.extend(taken);
        if before.is_some_and(|before| held < before) {
            break;
        }
        before = Some(held);
    }
    promoted
}
