use std::collections::{BinaryHeap, HashSet};

use log::{debug, warn};
use rayon::prelude::*;

use crate::embedding::knn::KnnGraph;
use crate::quality::Quality;
use crate::refine::dissolve::{
    DissolveLedger, DissolveSettings, POOL_VIEWS, PoolRun, PoolSearcher, PoolView, Pot, floor_for,
    neighbours_for,
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

fn rungs_of(
    pot: &Pot,
    pool: &HashSet<usize>,
    view: PoolView,
    settings: DissolveSettings,
    first: Option<&Built>,
    ledger: &mut DissolveLedger,
    search: &dyn PoolSearcher,
) -> Option<(Vec<Vec<usize>>, Built)> {
    let built = match first.and_then(|first| reuse(pool, first)) {
        Some(built) => built,
        None => match search.neighbours(pool, settings.n_neighbours, view) {
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
        ledger.search.linkage += found.len();
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
        let truncated;
        let graph = match width == built.knn.indices.ncols() {
            true => &built.knn,
            false => {
                truncated = built.knn.truncate(width);
                &truncated
            }
        };
        let results = match search.partition(graph, &built.order, round) {
            Ok(results) => results,
            Err(error) => {
                warn!("Could not partition the pool: {error}");
                break;
            }
        };
        ledger.search.rounds += 1;
        for result in results {
            ledger.search.noise += result.outliers.len();
            candidates.extend(result.cluster_map.into_values().map(sorted));
        }
    }
    Some((candidates, built))
}

fn propose(
    pot: &Pot,
    pool: &HashSet<usize>,
    settings: DissolveSettings,
    oracle: &[Vec<usize>],
    first: &mut Vec<Built>,
    ledger: &mut DissolveLedger,
    search: &dyn PoolSearcher,
) -> Vec<Vec<usize>> {
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
    if let [(_, combined), (PoolView::Composition, composition)] = per_view.as_slice() {
        ledger.search.composition += composition
            .iter()
            .filter(|group| combined.binary_search(group).is_err())
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
type Candidate = Ranked<Quality>;

fn heap(pot: &Pot, candidates: Vec<Vec<usize>>) -> BinaryHeap<Candidate> {
    candidates
        .into_par_iter()
        .map(|contigs| {
            let quality = pot.quality_of(&contigs);
            Ranked {
                worth: pot.worth_of(quality),
                contigs,
                extra: quality,
            }
        })
        .collect::<Vec<_>>()
        .into()
}

struct Taken {
    contigs: Vec<usize>,
    worth: f64,
}

struct Swept {
    taken: Vec<Taken>,
    refused: BinaryHeap<Candidate>,
    refusals: Vec<Refusal>,
    consumed: usize,
}

fn sweep(
    pot: &Pot,
    mut held: BinaryHeap<Candidate>,
    pool: &HashSet<usize>,
    claimed: &mut HashSet<usize>,
    bar: Rung,
    watch: Watch<'_, '_>,
) -> Swept {
    let mut swept = Swept {
        taken: Vec::new(),
        refused: BinaryHeap::new(),
        refusals: Vec::new(),
        consumed: 0,
    };
    while let Some(entry) = held.pop() {
        // A proposal the winners already emptied is not the proposal that was scored, so what is
        // left of it is put back through the same bar rather than trusted on the rank it earned.
        let left = remaining(&entry.contigs, claimed);
        if left.len() < 2 {
            swept.consumed += 1;
            watch.row(entry.worth, "consumed", &entry.contigs, pot, None);
            continue;
        }
        let quality = match left.len() == entry.contigs.len() {
            true => entry.extra,
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
                    swept.taken.push(Taken {
                        contigs: left,
                        worth,
                    });
                    continue;
                }
            }
            other => {
                watch.row(worth, other.label(), &left, pot, Some(quality));
                Refusal::Judged(other)
            }
        };
        // Ranking the remainder by its own worth here was measured and lost CAMI III bins.
        swept.refused.push(Ranked {
            worth: entry.worth,
            contigs: left,
            extra: quality,
        });
        swept.refusals.push(refusal);
    }
    swept
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

fn tally(refusals: &[Refusal], ledger: &mut DissolveLedger) {
    let refused = &mut ledger.refused;
    for refusal in refusals {
        match refusal {
            Refusal::Judged(Verdict::TooSmall) => refused.small += 1,
            Refusal::Judged(Verdict::Incomplete) => refused.incomplete += 1,
            Refusal::Judged(Verdict::Contaminated) => refused.contaminated += 1,
            Refusal::Judged(_) | Refusal::Worse => refused.worse += 1,
            Refusal::Carved => refused.carved += 1,
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
) -> (Vec<Taken>, Vec<Deferred>) {
    let mut promoted = Vec::new();
    let mut deferred = Vec::new();
    let mut claimed = HashSet::new();
    let mut held = heap(pot, candidates);
    let mut refusals = Vec::new();

    for at in 0..RUNGS {
        run.ledger.search.rung = run.ledger.search.rung.max(at);
        let bar = run.settings.bars.at(at);
        let watch = Watch {
            report: run.report,
            pass,
            rung: at,
        };
        let swept = sweep(pot, held, pool, &mut claimed, bar, watch);
        run.ledger.refused.consumed += swept.consumed;
        let empty = swept.taken.is_empty();
        match at > 0 {
            true => deferred.extend(swept.taken.into_iter().map(|taken| Deferred {
                contigs: taken.contigs,
                rung: at,
            })),
            false => promoted.extend(swept.taken),
        }
        held = swept.refused;
        refusals = swept.refusals;
        // Walking adopts at looser bars what the strict bar left, which a fragmented assembly
        // needs and a near-complete one loses bins to.
        if run.ledger.finished().mostly() && !empty {
            break;
        }
    }
    tally(&refusals, run.ledger);
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
        let swept = sweep(pot, heap(pot, candidates), pool, &mut claimed, bar, watch);
        run.ledger.refused.consumed += swept.consumed;
        run.ledger.claims.drained += swept.taken.len();
        promoted.extend(swept.taken.into_iter().map(|taken| taken.contigs));
        tally(&swept.refusals, run.ledger);
    }
    promoted
}

fn leave(pool: &mut HashSet<usize>, taken: &[Vec<usize>]) {
    for contig in taken.iter().flatten() {
        pool.remove(contig);
    }
}

// Every round searches the same pool, so the bar is asked which proposal to keep rather than
// which came first, and each pass then re-embeds what the pass before it left.
pub fn ranked(
    pot: &Pot,
    pool: &mut HashSet<usize>,
    run: &mut PoolRun<'_, '_>,
    oracle: &[Vec<usize>],
    search: &dyn PoolSearcher,
) -> Vec<Vec<usize>> {
    let settings = run.settings;
    let progress = crate::progress::hold(crate::progress::counted(
        crate::progress::Stage::RescuingUnbinned,
        settings.passes.max(1) as u64,
    ));
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
        run.ledger.search.proposed += candidates.len();

        let (taken, held) = {
            let _timer = crate::timing::scope("claim");
            claim(pot, candidates, &*pool, run, pass)
        };
        run.ledger.claims.deferred += held.len();
        deferred.extend(held);
        // The deferred keep their contigs in the pool, so a pass the strict bar takes nothing
        // from hands the next one the pool it just searched and would find the same clusters.
        if taken.is_empty() {
            break;
        }
        run.ledger.search.passes += 1;
        progress.inc(1);
        // How many passes a pool is worth differs per assembly, and the bins a pass finds are
        // worth less than the last one's long before it finds none, which no tier can see.
        let mean_worth = taken.iter().map(|taken| taken.worth).sum::<f64>() / taken.len() as f64;
        let taken = taken
            .into_iter()
            .map(|taken| taken.contigs)
            .collect::<Vec<_>>();
        leave(pool, &taken);
        promoted.extend(taken);
        if before.is_some_and(|before| mean_worth < before) {
            break;
        }
        before = Some(mean_worth);
    }
    crate::progress::done(&progress);
    if !deferred.is_empty() {
        let taken = drain(pot, deferred, pool, run, settings.passes.max(1));
        leave(pool, &taken);
        promoted.extend(taken);
    }
    promoted
}
