use std::collections::{BinaryHeap, HashSet};

use anyhow::Result;
use log::{debug, warn};
use rayon::prelude::*;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::knn::KnnGraph;
use crate::quality::Quality;
use crate::refine::dissolve::{
    DissolveLedger, DissolveSettings, POOL_VIEWS, PoolRun, PoolSearch, PoolView, Pot, RoundParams,
    floor_for, neighbours_for,
};
use crate::refine::pool_report::PoolReport;
use crate::refine::ranking::{Ranked, remaining, remaining_in, sorted};
use crate::refine::rung::{RUNGS, Rung, Verdict};

struct Built {
    knn: KnnGraph,
    order: Vec<usize>,
}

// A later pass searches a strict subset of the first one's pool, so its neighbours are already
// in that build and only the contigs that left have to be taken out of the rows.
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
    debug!(
        "Reused the pool graph for {} of {} contigs",
        keep.len(),
        first.order.len()
    );
    let order = keep.iter().map(|position| first.order[*position]).collect();
    Some(Built { knn, order })
}

fn rungs_of<N, P>(
    pot: &Pot,
    pool: &HashSet<usize>,
    view: PoolView,
    settings: DissolveSettings,
    first: Option<&Built>,
    ledger: &mut DissolveLedger,
    search: &PoolSearch<N, P>,
) -> Option<(Vec<Vec<usize>>, Built)>
where
    N: Fn(&HashSet<usize>, usize, PoolView) -> Result<(KnnGraph, Vec<usize>)>,
    P: Fn(&KnnGraph, &[usize], RoundParams) -> Result<Vec<Partitioning>>,
{
    let built = match first.and_then(|first| reuse(pool, first)) {
        Some(built) => built,
        None => match (search.neighbours)(pool, settings.n_neighbours, view) {
            Ok((knn, order)) => {
                debug!(
                    "Built the pool graph for {} contigs from scratch",
                    pool.len()
                );
                Built { knn, order }
            }
            Err(error) => {
                warn!("Could not re-embed the pool: {error}");
                return None;
            }
        },
    };
    let mut candidates = Vec::new();
    {
        let _timer = crate::timing::scope("linkage");
        let found = crate::refine::linkage::candidates(
            &built.knn,
            &built.order,
            |contig| pot.length(contig),
            floor_for(settings),
            settings.max_bin_size,
        );
        ledger.proposed_linkage += found.len();
        candidates.extend(found);
    }
    let mut last = 0;
    for round in 0..settings.rounds.max(1) {
        let round = neighbours_for(settings, round);
        // Every round reads the same build and differs only in how many neighbours it takes, so
        // two rounds the build cannot tell apart would run the same partition twice.
        let width = round.n_neighbours.min(built.knn.indices.ncols());
        if width == last {
            continue;
        }
        last = width;
        let truncated = built.knn.truncate(width);
        let results = match (search.partition)(&truncated, &built.order, round) {
            Ok(results) => results,
            Err(error) => {
                warn!("Could not partition the pool: {error}");
                break;
            }
        };
        ledger.rounds += 1;
        for result in results {
            ledger.noise += result.outliers.len();
            candidates.extend(result.cluster_map.into_values().map(sorted));
        }
    }
    Some((candidates, built))
}

fn propose<N, P>(
    pot: &Pot,
    pool: &HashSet<usize>,
    settings: DissolveSettings,
    oracle: &[Vec<usize>],
    first: &mut Vec<Built>,
    ledger: &mut DissolveLedger,
    search: &PoolSearch<N, P>,
) -> Vec<Vec<usize>>
where
    N: Fn(&HashSet<usize>, usize, PoolView) -> Result<(KnnGraph, Vec<usize>)>,
    P: Fn(&KnnGraph, &[usize], RoundParams) -> Result<Vec<Partitioning>>,
{
    let mut candidates = oracle
        .iter()
        .map(|group| remaining_in(group, pool))
        .filter(|group| group.len() >= 2)
        .collect::<Vec<_>>();
    let mut per_view = Vec::new();
    for (index, view) in POOL_VIEWS.iter().enumerate() {
        let Some((mut found, built)) =
            rungs_of(pot, pool, *view, settings, first.get(index), ledger, search)
        else {
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

// Only a candidate's last refusal is counted, so one refused at every rung counts once.
#[derive(Clone, Copy)]
enum Refusal {
    Judged(Verdict),
    Worse,
    Carved,
}

// Each entry carries the quality of the contigs it holds, scored here in parallel, so the
// serial sweep only scores again what an earlier claim took contigs from.
#[derive(Clone, Copy)]
struct Scored {
    refusal: Refusal,
    quality: Quality,
}

fn heap(pot: &Pot, candidates: Vec<Vec<usize>>) -> BinaryHeap<Ranked<Scored>> {
    candidates
        .into_par_iter()
        .map(|contigs| {
            let quality = pot.quality_of(&contigs);
            Ranked {
                worth: pot.worth_of(quality),
                contigs,
                extra: Scored {
                    refusal: Refusal::Worse,
                    quality,
                },
            }
        })
        .collect::<Vec<_>>()
        .into()
}

struct Swept {
    taken: Vec<Vec<usize>>,
    refused: BinaryHeap<Ranked<Scored>>,
    consumed: usize,
}

fn sweep(
    pot: &Pot,
    mut held: BinaryHeap<Ranked<Scored>>,
    pool: &HashSet<usize>,
    claimed: &mut HashSet<usize>,
    bar: Rung,
    watch: Watch<'_, '_>,
) -> Swept {
    let mut taken = Vec::new();
    let mut refused = BinaryHeap::new();
    let mut consumed = 0;
    while let Some(entry) = held.pop() {
        let left = remaining(&entry.contigs, claimed);
        if left.len() < 2 {
            consumed += 1;
            watch.row(entry.worth, "consumed", &entry.contigs, pot, None);
            continue;
        }
        let quality = match left.len() == entry.contigs.len() {
            true => entry.extra.quality,
            false => pot.quality_of(&left),
        };
        let worth = pot.worth_of(quality);
        let refusal = match pot.verdict(&left, quality, bar) {
            Verdict::Adopt => {
                let standing = pot.standing(&left, pool, claimed);
                if !pot.conserves(&left, worth, &standing) {
                    watch.row(worth, "worse", &left, pot, Some(quality));
                    Refusal::Worse
                } else if !pot.unifies(quality, &standing) {
                    watch.row(worth, "carves", &left, pot, Some(quality));
                    Refusal::Carved
                } else {
                    claimed.extend(left.iter().copied());
                    watch.row(worth, Verdict::Adopt.label(), &left, pot, Some(quality));
                    taken.push(left);
                    continue;
                }
            }
            other => {
                watch.row(worth, other.label(), &left, pot, Some(quality));
                Refusal::Judged(other)
            }
        };
        // Ranking the remainder by its own worth here was measured and lost CAMI III bins.
        refused.push(Ranked {
            worth: entry.worth,
            contigs: left,
            extra: Scored { refusal, quality },
        });
    }
    Swept {
        taken,
        refused,
        consumed,
    }
}

#[derive(Clone, Copy)]
struct Watch<'a, 'n> {
    report: Option<&'a PoolReport<'n>>,
    pass: usize,
    rung: usize,
}

impl Watch<'_, '_> {
    fn row(
        &self,
        worth: f64,
        verdict: &str,
        contigs: &[usize],
        pot: &Pot,
        quality: Option<Quality>,
    ) {
        let Some(report) = self.report else {
            return;
        };
        report.row(crate::refine::pool_report::Row {
            pass: self.pass,
            rung: self.rung,
            worth,
            bp: pot.bases(contigs),
            quality: quality.unwrap_or_else(|| pot.quality_of(contigs)),
            verdict,
            contigs,
            origins: &pot.origins(contigs),
        });
    }
}

fn tally(refused: &BinaryHeap<Ranked<Scored>>, ledger: &mut DissolveLedger) {
    for entry in refused {
        match entry.extra.refusal {
            Refusal::Judged(Verdict::TooSmall) => ledger.refused_small += 1,
            Refusal::Judged(Verdict::Incomplete) => ledger.refused_incomplete += 1,
            Refusal::Judged(Verdict::Contaminated) => ledger.refused_contaminated += 1,
            Refusal::Judged(_) | Refusal::Worse => ledger.refused_worse += 1,
            Refusal::Carved => ledger.refused_carved += 1,
        }
    }
}

struct Deferred {
    contigs: Vec<usize>,
    rung: usize,
}

fn claim(
    pot: &Pot,
    candidates: Vec<Vec<usize>>,
    pool: &HashSet<usize>,
    run: &mut PoolRun<'_, '_>,
    pass: usize,
) -> (Vec<Vec<usize>>, Vec<Deferred>) {
    let mut promoted = Vec::new();
    let mut deferred = Vec::new();
    let mut claimed = HashSet::new();
    let mut held = heap(pot, candidates);

    for at in 0..RUNGS {
        run.ledger.rung = run.ledger.rung.max(at);
        let bar = run.settings.bars.at(at);
        let watch = Watch {
            report: run.report,
            pass,
            rung: at,
        };
        let swept = sweep(pot, held, pool, &mut claimed, bar, watch);
        let Swept {
            taken,
            refused,
            consumed,
        } = swept;
        run.ledger.refused_consumed += consumed;
        let empty = taken.is_empty();
        match at > 0 {
            true => deferred.extend(
                taken
                    .into_iter()
                    .map(|contigs| Deferred { contigs, rung: at }),
            ),
            false => promoted.extend(taken),
        }
        held = refused;
        // Walking adopts at looser bars what the strict bar left, which a fragmented assembly
        // needs and a near-complete one loses bins to.
        if run.ledger.finished().mostly() && !empty {
            break;
        }
    }
    tally(&held, run.ledger);
    (promoted, deferred)
}

// Settled against the pool the passes left, not the one they started on, so a loose rung only
// gets what the strict bar had every pass to want and did not take.
fn drain(
    pot: &Pot,
    deferred: Vec<Deferred>,
    pool: &HashSet<usize>,
    run: &mut PoolRun<'_, '_>,
    pass: usize,
) -> Vec<Vec<usize>> {
    let mut by_rung = vec![Vec::new(); RUNGS];
    for entry in deferred {
        let left = remaining_in(&entry.contigs, pool);
        if left.len() >= 2 {
            by_rung[entry.rung].push(left);
        }
    }
    let mut promoted = Vec::new();
    let mut claimed = HashSet::new();
    for (at, slot) in by_rung.iter_mut().enumerate().skip(1) {
        let mut candidates = std::mem::take(slot);
        dedupe(&mut candidates);
        if candidates.is_empty() {
            continue;
        }
        let bar = run.settings.bars.at(at);
        let watch = Watch {
            report: run.report,
            pass,
            rung: at,
        };
        let Swept {
            taken,
            refused,
            consumed,
        } = sweep(pot, heap(pot, candidates), pool, &mut claimed, bar, watch);
        run.ledger.refused_consumed += consumed;
        run.ledger.drained += taken.len();
        promoted.extend(taken);
        tally(&refused, run.ledger);
    }
    promoted
}

// Every round searches the same pool, so the bar is asked which proposal to keep rather than
// which came first, and each pass then re-embeds what the pass before it left.
pub fn ranked<N, P>(
    pot: &Pot,
    pool: &mut HashSet<usize>,
    run: &mut PoolRun<'_, '_>,
    oracle: &[Vec<usize>],
    search: &PoolSearch<N, P>,
) -> Vec<Vec<usize>>
where
    N: Fn(&HashSet<usize>, usize, PoolView) -> Result<(KnnGraph, Vec<usize>)>,
    P: Fn(&KnnGraph, &[usize], RoundParams) -> Result<Vec<Partitioning>>,
{
    let settings = run.settings;
    let progress = crate::progress::counted(
        crate::progress::Stage::RescuingUnbinned,
        settings.passes.max(1) as u64,
    );
    let mut promoted = Vec::new();
    let mut deferred = Vec::new();
    let mut before: Option<f64> = None;
    let mut first = Vec::new();
    for pass in 0..settings.passes.max(1) {
        progress.set_message(format!("{} in the pool", pool.len()));
        if pool.len() < crate::refine::dissolve::MIN_RESCUE_CONTIGS {
            break;
        }
        let mut candidates = propose(pot, pool, settings, oracle, &mut first, run.ledger, search);
        dedupe(&mut candidates);
        run.ledger.proposed += candidates.len();

        let (taken, held) = {
            let _timer = crate::timing::scope("claim");
            claim(pot, candidates, &*pool, run, pass)
        };
        run.ledger.deferred += held.len();
        deferred.extend(held);
        // The deferred keep their contigs in the pool, so a pass the strict bar takes nothing
        // from hands the next one the pool it just searched and would find the same clusters.
        if taken.is_empty() {
            break;
        }
        for contigs in &taken {
            for contig in contigs {
                pool.remove(contig);
            }
        }
        run.ledger.passes += 1;
        progress.inc(1);
        // How many passes a pool is worth differs per assembly, and the bins a pass finds are
        // worth less than the last one's long before it finds none, which no tier can see.
        let mean_worth =
            taken.iter().map(|contigs| pot.worth(contigs)).sum::<f64>() / taken.len() as f64;
        promoted.extend(taken);
        if before.is_some_and(|before| mean_worth < before) {
            break;
        }
        before = Some(mean_worth);
    }
    progress.finish_and_clear();
    if !deferred.is_empty() {
        let pass = settings.passes.max(1);
        let taken = drain(pot, deferred, pool, run, pass);
        for contigs in &taken {
            for contig in contigs {
                pool.remove(contig);
            }
        }
        promoted.extend(taken);
    }
    promoted
}
