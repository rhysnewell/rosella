use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};

use anyhow::Result;
use log::warn;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::knn::KnnGraph;
use crate::refine::dissolve::{
    DissolveLedger, DissolveSettings, POOL_VIEWS, PoolView, Pot, RoundParams, neighbours_for,
};
use crate::refine::rung::{RUNGS, Rung, Verdict};

pub struct Built {
    knn: KnnGraph,
    order: Vec<usize>,
}

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

/// A later pass searches a strict subset of the first one's pool, so its neighbours are already
/// in that build and only the contigs that left have to be taken out of the rows.
fn reuse(pool: &HashSet<usize>, first: &Built) -> Option<Built> {
    let keep = first
        .order
        .iter()
        .enumerate()
        .filter(|(_, contig)| pool.contains(contig))
        .map(|(position, _)| position)
        .collect::<Vec<_>>();
    if keep.len() == first.order.len() {
        return None;
    }
    let knn = first.knn.induced(&keep)?;
    let order = keep.iter().map(|position| first.order[*position]).collect();
    Some(Built { knn, order })
}

fn rungs_of(
    pool: &HashSet<usize>,
    view: PoolView,
    settings: DissolveSettings,
    first: Option<&Built>,
    ledger: &mut DissolveLedger,
    neighbours: &impl Fn(&HashSet<usize>, usize, PoolView) -> Result<(KnnGraph, Vec<usize>)>,
    partition: &impl Fn(&KnnGraph, &[usize], RoundParams) -> Result<Vec<Partitioning>>,
) -> Option<(Vec<Vec<usize>>, Built)> {
    let built = match first.filter(|_| settings.reuse).and_then(|first| reuse(pool, first)) {
        Some(built) => built,
        None => match neighbours(pool, settings.n_neighbours, view) {
            Ok((knn, order)) => Built { knn, order },
            Err(error) => {
                warn!("Could not re-embed the pool: {error}");
                return None;
            }
        },
    };
    let mut candidates = Vec::new();
    let mut last = 0;
    for round in 0..settings.rounds.max(1) {
        let round = neighbours_for(settings, round);
        // Every round reads the same build and differs only in how many neighbours it takes, so
        // two rounds the build cannot tell apart would run the same partition twice.
        let truncated = built.knn.truncate(round.n_neighbours);
        let width = truncated.indices.ncols();
        if width == last {
            continue;
        }
        last = width;
        let results = match partition(&truncated, &built.order, round) {
            Ok(results) => results,
            Err(error) => {
                warn!("Could not partition the pool: {error}");
                break;
            }
        };
        ledger.rounds += 1;
        for result in results {
            ledger.noise = result.outliers.len();
            candidates.extend(result.cluster_map.into_values().map(sorted));
        }
    }
    Some((candidates, built))
}

fn propose(
    pool: &HashSet<usize>,
    settings: DissolveSettings,
    oracle: &[Vec<usize>],
    first: &mut Vec<Built>,
    ledger: &mut DissolveLedger,
    neighbours: &impl Fn(&HashSet<usize>, usize, PoolView) -> Result<(KnnGraph, Vec<usize>)>,
    partition: &impl Fn(&KnnGraph, &[usize], RoundParams) -> Result<Vec<Partitioning>>,
) -> Vec<Vec<usize>> {
    let mut candidates = oracle
        .iter()
        .map(|group| remaining_in(group, pool))
        .filter(|group| group.len() >= 2)
        .collect::<Vec<_>>();
    let mut per_view = Vec::new();
    for (index, view) in POOL_VIEWS.iter().enumerate() {
        let Some((mut found, built)) = rungs_of(
            pool,
            *view,
            settings,
            first.get(index),
            ledger,
            neighbours,
            partition,
        ) else {
            continue;
        };
        if first.len() == index {
            first.push(built);
        }
        dedupe(&mut found);
        per_view.push((*view, found));
    }
    if let [(_, first), (PoolView::Composition, second)] = per_view.as_slice() {
        ledger.proposed_composition += second
            .iter()
            .filter(|group| first.binary_search(group).is_err())
            .count();
    }
    for (_, found) in per_view {
        candidates.extend(found);
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
    neighbours: impl Fn(&HashSet<usize>, usize, PoolView) -> Result<(KnnGraph, Vec<usize>)>,
    partition: impl Fn(&KnnGraph, &[usize], RoundParams) -> Result<Vec<Partitioning>>,
) -> Vec<Vec<usize>> {
    let mut promoted = Vec::new();
    let mut before: Option<f64> = None;
    let mut first = Vec::new();
    for _ in 0..settings.passes.max(1) {
        if pool.len() < settings.min_contigs {
            break;
        }
        let mut candidates =
            propose(pool, settings, oracle, &mut first, ledger, &neighbours, &partition);
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
