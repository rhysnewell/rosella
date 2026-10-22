use std::collections::{BTreeMap, HashSet};

use anyhow::Result;
use log::{debug, warn};

use crate::clustering::clusterer::{Partitioning, placed_once};
use crate::embedding::features::ContigFeatures;
use crate::embedding::knn::KnnGraph;
use crate::quality::Scorer;
use crate::refine::pool_report::{PoolReport, Row};
use crate::refine::ranking::{remaining, sorted};
use crate::refine::rung::{Bars, Verdict, verdict};
use crate::refine::select::ranked;

mod ledger;
mod pot;

pub use ledger::{DissolveLedger, Mass};
pub use pot::{Pot, Standing};

const MIN_NEIGHBOURS: usize = 2;

// Every round searches the whole pool at half the neighbours of the one before it, so a genome
// the dense graph buries can still form its own community in a sparser one.
#[derive(Debug, Clone, Copy)]
pub struct RoundParams {
    pub n_neighbours: usize,
    pub ladder: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
// At one sample the coverage half of the distance is a tie on depth alone, so the two views
// hold different genomes and the bar is asked which grouping to keep rather than which metric.
pub enum PoolView {
    Combined,
    Composition,
}

pub const POOL_VIEWS: [PoolView; 2] = [PoolView::Combined, PoolView::Composition];

// At weight zero the combined distance is the composition one, so a second view repeats it.
fn pool_views(features: &ContigFeatures) -> &'static [PoolView] {
    match features.distance_settings().is_composition_only() {
        true => &POOL_VIEWS[..1],
        false => &POOL_VIEWS,
    }
}

// The fuzzy set needs two neighbours, and a subset of three is the smallest that has them.
pub(crate) const MIN_RESCUE_CONTIGS: usize = 3;

#[derive(Debug, Clone, Copy)]
pub struct DissolveSettings {
    pub bars: Bars,
    pub rounds: usize,
    pub passes: usize,
    pub n_neighbours: usize,
    pub max_bin_size: usize,
}

pub struct PoolInputs<'a, 'n> {
    pub features: &'a ContigFeatures<'a>,
    pub quality: &'a dyn Scorer,
    pub settings: DissolveSettings,
    pub oracle: &'a [Vec<usize>],
    pub report: Option<&'a PoolReport<'n>>,
}

pub trait PoolSearcher {
    fn neighbours(
        &self,
        pool: &HashSet<usize>,
        n_neighbours: usize,
        view: PoolView,
    ) -> Result<(KnnGraph, Vec<usize>)>;

    fn partition(
        &self,
        knn: &KnnGraph,
        order: &[usize],
        round: RoundParams,
    ) -> Result<Vec<Partitioning>>;
}

pub struct PoolSearch<N, P> {
    pub neighbours: N,
    pub partition: P,
}

impl<N, P> PoolSearch<N, P>
where
    N: Fn(&HashSet<usize>, usize, PoolView) -> Result<(KnnGraph, Vec<usize>)>,
    P: Fn(&KnnGraph, &[usize], RoundParams) -> Result<Vec<Partitioning>>,
{
    pub fn new(neighbours: N, partition: P) -> Self {
        Self {
            neighbours,
            partition,
        }
    }
}

impl<N, P> PoolSearcher for PoolSearch<N, P>
where
    N: Fn(&HashSet<usize>, usize, PoolView) -> Result<(KnnGraph, Vec<usize>)>,
    P: Fn(&KnnGraph, &[usize], RoundParams) -> Result<Vec<Partitioning>>,
{
    fn neighbours(
        &self,
        pool: &HashSet<usize>,
        n_neighbours: usize,
        view: PoolView,
    ) -> Result<(KnnGraph, Vec<usize>)> {
        (self.neighbours)(pool, n_neighbours, view)
    }

    fn partition(
        &self,
        knn: &KnnGraph,
        order: &[usize],
        round: RoundParams,
    ) -> Result<Vec<Partitioning>> {
        (self.partition)(knn, order, round)
    }
}

pub struct PoolRun<'a, 'n> {
    pub settings: DissolveSettings,
    pub views: &'static [PoolView],
    pub ledger: &'a mut DissolveLedger,
    pub report: Option<&'a PoolReport<'n>>,
}

// A bin the pool cannot be expected to improve is held out, because the pool re-partitions what
// it is handed and on a strain-heavy assembly that bisects whole genomes into two half bins.
fn dissolving(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    bins: &BTreeMap<usize, Vec<usize>>,
    settings: DissolveSettings,
    report: Option<&PoolReport<'_>>,
    ledger: &mut DissolveLedger,
) -> Vec<(usize, Vec<usize>)> {
    let top = settings.bars.piece();
    let bar = settings.bars.at(0);
    let mut dissolving = Vec::new();
    for (bin_id, contigs) in bins.iter() {
        let size = features.bin_size(contigs);
        let scored = quality.score(contigs);
        let held = verdict(size, scored, bar) == Verdict::Adopt;
        if let Some(report) = report {
            report.row(Row {
                pass: 0,
                rung: 0,
                worth: scored.score(settings.bars.worth),
                bp: size,
                quality: scored,
                verdict: if held { "held" } else { "dissolved" },
                contigs,
                origins: &[(*bin_id, size)],
            });
        }
        let handed = &mut ledger.handed;
        if held {
            handed.held_back += 1;
            continue;
        }
        if size < top {
            handed.small += 1;
        } else {
            handed.clean += 1;
        }
        handed.bp += size;
        dissolving.push((*bin_id, contigs.clone()));
    }
    dissolving
}

// The probe asks whether the bar takes the right grouping when it is handed one, so what
// matters is how many of the offered groups came back out, not how many were proposed.
fn report_oracle(
    features: &ContigFeatures,
    handed: &HashSet<usize>,
    oracle: &[Vec<usize>],
    promoted: &[Vec<usize>],
) {
    let taken = promoted.iter().collect::<HashSet<_>>();
    let offered = oracle
        .iter()
        .map(|group| {
            group
                .iter()
                .copied()
                .filter(|contig| handed.contains(contig))
                .collect::<Vec<_>>()
        })
        .filter(|group| group.len() >= 2)
        .collect::<Vec<_>>();
    let whole = offered.iter().filter(|group| taken.contains(group)).count();
    let bp = offered
        .iter()
        .filter(|group| taken.contains(group))
        .map(|group| features.bin_size(group))
        .sum::<usize>();
    debug!(
        "Oracle groups: {} offered, {} taken whole ({} bp)",
        offered.len(),
        whole,
        bp
    );
}

pub fn neighbours_for(settings: DissolveSettings, round: usize) -> RoundParams {
    RoundParams {
        n_neighbours: settings
            .n_neighbours
            .checked_shr(round as u32)
            .unwrap_or(0)
            .max(MIN_NEIGHBOURS),
        ladder: true,
    }
}

// Every bin goes back in the pot with the unbinned and is searched again without the bins that
// already left, which is the one thing re-cutting inside a bin cannot do.
pub fn dissolve(
    inputs: PoolInputs<'_, '_>,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    unbinned: &mut Vec<usize>,
    search: &dyn PoolSearcher,
) -> DissolveLedger {
    let PoolInputs {
        features,
        quality,
        settings,
        oracle,
        report,
    } = inputs;
    let mut ledger = DissolveLedger::default();
    let dissolved = dissolving(features, quality, bins, settings, report, &mut ledger);

    let mut pool = unbinned.iter().copied().collect::<HashSet<_>>();
    for (_, contigs) in &dissolved {
        pool.extend(contigs.iter().copied());
    }
    ledger.pool = Mass::of(features, &pool);

    if pool.len() < MIN_RESCUE_CONTIGS {
        return ledger;
    }
    let handed = pool.clone();
    let pot = Pot::new(
        features,
        quality,
        settings.bars.worth,
        settings.bars.piece(),
        &dissolved,
    );
    let mut run = PoolRun {
        settings,
        views: pool_views(features),
        ledger: &mut ledger,
        report,
    };
    let promoted = ranked(&pot, &mut pool, &mut run, oracle, search);
    if !oracle.is_empty() {
        report_oracle(features, &handed, oracle, &promoted);
    }
    let judge = crate::refine::restore::Judge {
        features,
        quality,
        reported: settings.bars.reported(),
        accept: settings.bars.at(0),
    };
    let held = crate::refine::restore::restore(&judge, settings.bars.worth, &dissolved, promoted);
    ledger.claims.restored = held.bins;
    pool.extend(held.released);
    let mut promoted = held.promoted;
    if promoted.is_empty() {
        return ledger;
    }
    let folded = crate::refine::foldback::fold_back(
        features,
        quality,
        settings.bars.worth,
        settings.bars.reported(),
        &dissolved,
        &mut promoted,
    );
    ledger.placed.folded = Mass::of(features, &folded);
    for contig in &folded {
        pool.remove(contig);
    }
    // Bin ids follow this order, so it rests on each bin's first contig.
    promoted.sort_unstable_by_key(|contigs| contigs.first().copied().unwrap_or(usize::MAX));

    if let Err(error) = placed_once(promoted.iter().flatten().copied(), &handed) {
        warn!("Pool: {error} (bins unchanged)");
        return ledger;
    }

    let next_bin_id = bins.keys().max().map_or(0, |id| id + 1);
    let mut claimed = HashSet::new();
    ledger.claims.promoted = promoted.len();
    for (offset, contigs) in promoted.into_iter().enumerate() {
        ledger.placed.adopted.add(Mass::of(features, &contigs));
        claimed.extend(contigs.iter().copied());
        bins.insert(next_bin_id + offset, contigs);
    }

    for (bin_id, contigs) in dissolved {
        let kept = remaining(&contigs, &claimed);
        if kept.is_empty() {
            ledger.placed.emptied += 1;
            bins.remove(&bin_id);
            continue;
        }
        ledger.placed.returned.add(Mass::of(features, &kept));
        for contig in &kept {
            pool.remove(contig);
        }
        bins.insert(bin_id, kept);
    }

    ledger.placed.left = Mass::of(features, &pool);
    *unbinned = sorted(pool);
    ledger
}
