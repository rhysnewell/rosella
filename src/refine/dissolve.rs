use std::collections::{BTreeMap, HashMap, HashSet};

use anyhow::Result;
use log::{debug, warn};

use crate::clustering::clusterer::{Partitioning, placed_once};
use crate::embedding::features::ContigFeatures;
use crate::embedding::knn::KnnGraph;
use crate::quality::Scorer;
use crate::refine::finished::Finished;
use crate::refine::pool_report::PoolReport;
use crate::refine::rung::{Bars, Rung, Verdict, judge};
use crate::refine::select::{ranked, remaining, remaining_in, sorted};

const MIN_NEIGHBOURS: usize = 2;

/// Every round searches the whole pool at half the neighbours of the one before it, so a genome
/// the dense graph buries can still form its own community in a sparser one.
#[derive(Debug, Clone, Copy)]
pub struct RoundParams {
    pub n_neighbours: usize,
    pub ladder: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// At one sample the coverage half of the distance is a tie on depth alone, so the two views
/// hold different genomes and the bar is asked which grouping to keep rather than which metric.
pub enum PoolView {
    Combined,
    Composition,
}

pub const POOL_VIEWS: [PoolView; 2] = [PoolView::Combined, PoolView::Composition];

/// What the pool keeps out of the pot. The bars are what a bin has to clear to be adopted,
/// which is a stricter question than whether re-partitioning it can do better.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Hold {
    Bars,
    Size,
    Tier,
}

impl Hold {
    pub fn holds(
        self,
        features: &ContigFeatures,
        quality: &dyn Scorer,
        contigs: &[usize],
        bars: Bars,
        rung: Rung,
    ) -> bool {
        match self {
            Self::Bars => judge(features, quality, contigs, rung) == Verdict::Adopt,
            _ if features.bin_size(contigs) < rung.floor => false,
            Self::Size => true,
            Self::Tier => quality.score(contigs).contamination <= bars.tier(),
        }
    }
}

/// How far a pass walks its rungs. Walking adopts at looser bars what the strict bar left,
/// which a fragmented assembly needs and a near-complete one loses bins to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RungWalk {
    Walk,
    Break,
}

#[derive(Debug, Clone, Copy)]
pub struct DissolveSettings {
    pub bars: Bars,
    pub hold: Hold,
    pub genome_floor: Option<usize>,
    pub min_contigs: usize,
    pub rounds: usize,
    pub passes: usize,
    pub n_neighbours: usize,
    pub max_bin_size: usize,
    pub reembed: bool,
    pub rung_walk: RungWalk,
    pub seed: u64,
}

pub struct PoolInputs<'a, 'n> {
    pub features: &'a ContigFeatures<'a>,
    pub quality: &'a dyn Scorer,
    pub settings: DissolveSettings,
    pub oracle: &'a [Vec<usize>],
    pub report: Option<&'a PoolReport<'n>>,
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

pub struct PoolRun<'a, 'n> {
    pub settings: DissolveSettings,
    pub top: usize,
    pub ledger: &'a mut DissolveLedger,
    pub report: Option<&'a PoolReport<'n>>,
}

/// What the pool took, what it refused and where the refusals went, in contigs and bases.
/// The bin count alone cannot say whether an arm has headroom left.
#[derive(Debug, Default, Clone, Copy)]
pub struct DissolveLedger {
    pub dissolved_small: usize,
    pub dissolved_clean: usize,
    pub dissolved_bp: usize,
    pub held_back: usize,
    pub pool_contigs: usize,
    pub pool_bp: usize,
    pub rounds: usize,
    pub passes: usize,
    pub rung: usize,
    pub proposed: usize,
    pub proposed_composition: usize,
    pub proposed_linkage: usize,
    pub refused_small: usize,
    pub refused_incomplete: usize,
    pub refused_contaminated: usize,
    pub refused_consumed: usize,
    pub refused_worse: usize,
    pub refused_carved: usize,
    pub noise: usize,
    pub promoted: usize,
    pub adopted_contigs: usize,
    pub adopted_bp: usize,
    pub returned_contigs: usize,
    pub returned_bp: usize,
    pub emptied: usize,
    pub restored: usize,
    pub deferred: usize,
    pub drained: usize,
    pub left_contigs: usize,
    pub left_bp: usize,
}

impl DissolveLedger {
    pub fn finished(&self) -> Finished {
        Finished::new(
            self.held_back,
            self.held_back + self.dissolved_small + self.dissolved_clean,
        )
    }
}

impl std::fmt::Display for DissolveLedger {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "held back {} bins already over the bars; dissolved {} small and {} \
             clean bins holding {} bp; pool {} contigs \
             {} bp; {} rounds over {} passes ending at rung {}; proposed {} clusters, refused {} small, {} \
             incomplete, {} contaminated, {} eaten by a rival, {} no better and {} carved out of \
             a bin with no duplication to explain them, {} noise; {} of the \
             proposals came only from composition, {} from the merge order; promoted {} \
             bins adopting {} contigs {} bp; returned {} contigs {} bp, emptied {} bins; left {} \
             contigs {} bp unbinned; restored {} bins the pool broke into nothing; held {} \
             loose candidates back and drained {} of them",
            self.held_back,
            self.dissolved_small,
            self.dissolved_clean,
            self.dissolved_bp,
            self.pool_contigs,
            self.pool_bp,
            self.rounds,
            self.passes,
            self.rung,
            self.proposed,
            self.refused_small,
            self.refused_incomplete,
            self.refused_contaminated,
            self.refused_consumed,
            self.refused_worse,
            self.refused_carved,
            self.noise,
            self.proposed_composition,
            self.proposed_linkage,
            self.promoted,
            self.adopted_contigs,
            self.adopted_bp,
            self.returned_contigs,
            self.returned_bp,
            self.emptied,
            self.left_contigs,
            self.left_bp,
            self.restored,
            self.deferred,
            self.drained
        )
    }
}

/// A piece smaller than the run's own genome scale is a shard of one, and relaxing that to let
/// more of the pool through cost more bins than it recovered.
pub(crate) fn floor_for(settings: DissolveSettings) -> usize {
    settings
        .genome_floor
        .unwrap_or(settings.bars.min_bin_size)
        .max(settings.bars.min_bin_size)
}

fn bases(features: &ContigFeatures, contigs: &HashSet<usize>) -> usize {
    contigs.iter().map(|contig| features.length(*contig)).sum()
}

/// A bin the pool cannot be expected to improve is held out, because the pool re-partitions what
/// it is handed and on a strain-heavy assembly that bisects whole genomes into two half bins.
fn dissolving(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    bins: &BTreeMap<usize, Vec<usize>>,
    top: usize,
    settings: DissolveSettings,
    report: Option<&crate::refine::pool_report::PoolReport<'_>>,
    ledger: &mut DissolveLedger,
) -> Vec<(usize, Vec<usize>)> {
    let bar = settings.bars.at(top, 0);
    let mut dissolving = Vec::new();
    for (bin_id, contigs) in bins.iter() {
        let held = settings
            .hold
            .holds(features, quality, contigs, settings.bars, bar);
        if let Some(report) = report {
            let scored = quality.score(contigs);
            let size = features.bin_size(contigs);
            report.row(crate::refine::pool_report::Row {
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
        if held {
            ledger.held_back += 1;
            continue;
        }
        if features.bin_size(contigs) < top {
            ledger.dissolved_small += 1;
        } else {
            ledger.dissolved_clean += 1;
        }
        ledger.dissolved_bp += features.bin_size(contigs);
        dissolving.push((*bin_id, contigs.clone()));
    }
    dissolving
}

pub struct Pot<'a> {
    features: &'a ContigFeatures<'a>,
    quality: &'a dyn Scorer,
    worth: f64,
    floor: usize,
    origin: HashMap<usize, usize>,
    members: HashMap<usize, Vec<usize>>,
}

impl<'a> Pot<'a> {
    pub fn new(
        features: &'a ContigFeatures<'a>,
        quality: &'a dyn Scorer,
        worth: f64,
        floor: usize,
        dissolved: &[(usize, Vec<usize>)],
    ) -> Self {
        Self {
            features,
            quality,
            worth,
            floor,
            origin: dissolved
                .iter()
                .flat_map(|(bin_id, contigs)| contigs.iter().map(|contig| (*contig, *bin_id)))
                .collect(),
            members: dissolved.iter().cloned().collect(),
        }
    }

    pub fn length(&self, contig: usize) -> usize {
        self.features.length(contig)
    }

    pub fn bases(&self, contigs: &[usize]) -> usize {
        self.features.bin_size(contigs)
    }

    pub fn quality_of(&self, contigs: &[usize]) -> crate::quality::Quality {
        self.quality.score(contigs)
    }

    pub fn worth(&self, contigs: &[usize]) -> f64 {
        self.quality.score(contigs).score(self.worth)
    }

    pub fn judge(&self, contigs: &[usize], rung: Rung) -> Verdict {
        judge(self.features, self.quality, contigs, rung)
    }

    pub fn origins(&self, contigs: &[usize]) -> Vec<(usize, usize)> {
        let mut taken: HashMap<usize, usize> = HashMap::new();
        for contig in contigs {
            if let Some(bin) = self.origin.get(contig) {
                *taken.entry(*bin).or_default() += self.features.length(*contig);
            }
        }
        let mut taken = taken.into_iter().collect::<Vec<_>>();
        taken.sort_unstable();
        taken
    }

    /// The remainder is where the loss sits. A candidate outscores the contigs it drains
    /// almost by construction, so what has to hold is that the bin left behind is no worse
    /// than the bin found, unless the candidate is itself at least that good.
    pub fn conserves(
        &self,
        contigs: &[usize],
        pool: &HashSet<usize>,
        claimed: &HashSet<usize>,
    ) -> bool {
        let candidate = self.worth(contigs);
        let taking = contigs.iter().copied().collect::<HashSet<_>>();
        self.origins(contigs).into_iter().all(|(bin, _)| {
            let Some(members) = self.members.get(&bin) else {
                return true;
            };
            let standing = remaining(&remaining_in(members, pool), claimed);
            let before = self.worth(&standing);
            candidate >= before || self.worth(&remaining(&standing, &taking)) >= before
        })
    }
}

/// A strain half reads complete and clean, so worth cannot tell a genome carved in two from an
/// organism pulled out of a bin holding two. Only the second leaves the duplication behind.
impl Pot<'_> {
    pub fn unifies(
        &self,
        contigs: &[usize],
        pool: &HashSet<usize>,
        claimed: &HashSet<usize>,
    ) -> bool {
        let taken = self.origins(contigs);
        let [(bin, bases)] = taken[..] else {
            return true;
        };
        let Some(members) = self.members.get(&bin) else {
            return true;
        };
        let standing = remaining(&remaining_in(members, pool), claimed);
        if self.bases(&standing).saturating_sub(bases) < self.floor {
            return true;
        }
        let doubled = self.quality_of(&standing).contamination;
        doubled > 0.0 && 2.0 * self.quality_of(contigs).contamination <= doubled
    }
}

/// The probe asks whether the bar takes the right grouping when it is handed one, so what
/// matters is how many of the offered groups came back out, not how many were proposed.
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
        "Oracle groups: offered {}, taken whole {} holding {} bp",
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

/// Every bin goes back in the pot with the unbinned and is searched again without the bins that
/// already left, which is the one thing re-cutting inside a bin cannot do.
pub fn dissolve<N, P>(
    inputs: PoolInputs<'_, '_>,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    unbinned: &mut Vec<usize>,
    search: PoolSearch<N, P>,
) -> DissolveLedger
where
    N: Fn(&HashSet<usize>, usize, PoolView) -> Result<(KnnGraph, Vec<usize>)>,
    P: Fn(&KnnGraph, &[usize], RoundParams) -> Result<Vec<Partitioning>>,
{
    let PoolInputs {
        features,
        quality,
        settings,
        oracle,
        report,
    } = inputs;
    let mut settings = settings;
    let mut ledger = DissolveLedger::default();
    let top = floor_for(settings);
    let dissolved = dissolving(features, quality, bins, top, settings, report, &mut ledger);
    if ledger.finished().mostly() {
        settings.rung_walk = RungWalk::Break;
    }

    let mut pool = unbinned.iter().copied().collect::<HashSet<_>>();
    for (_, contigs) in &dissolved {
        pool.extend(contigs.iter().copied());
    }
    ledger.pool_contigs = pool.len();
    ledger.pool_bp = bases(features, &pool);

    if pool.len() < settings.min_contigs {
        return ledger;
    }
    let handed = pool.clone();
    let pot = Pot::new(features, quality, settings.bars.worth, top, &dissolved);
    let mut run = PoolRun {
        settings,
        top,
        ledger: &mut ledger,
        report,
    };
    let promoted = ranked(&pot, &mut pool, &mut run, oracle, &search);
    if !oracle.is_empty() {
        report_oracle(features, &handed, oracle, &promoted);
    }
    let judge = crate::refine::restore::Judge {
        features,
        quality,
        reported: settings.bars.reported(top),
        accept: settings.bars.at(top, 0),
    };
    let held = crate::refine::restore::restore(&judge, settings.bars.worth, &dissolved, promoted);
    ledger.restored = held.bins;
    pool.extend(held.released);
    let mut promoted = held.promoted;
    if promoted.is_empty() {
        return ledger;
    }
    // Hash order would give the same partition different bin names on every run.
    promoted.sort_unstable_by_key(|contigs| contigs.first().copied().unwrap_or(usize::MAX));

    if let Err(error) = placed_once(promoted.iter().flatten().copied(), &handed) {
        warn!("The pool did not come back partitioned, leaving the bins alone: {error}");
        return ledger;
    }

    let next_bin_id = bins.keys().max().map_or(0, |id| id + 1);
    let mut claimed = HashSet::new();
    ledger.promoted = promoted.len();
    for (offset, contigs) in promoted.into_iter().enumerate() {
        ledger.adopted_contigs += contigs.len();
        ledger.adopted_bp += features.bin_size(&contigs);
        claimed.extend(contigs.iter().copied());
        bins.insert(next_bin_id + offset, contigs);
    }

    for (bin_id, contigs) in dissolved {
        let kept = remaining(&contigs, &claimed);
        if kept.is_empty() {
            ledger.emptied += 1;
            bins.remove(&bin_id);
            continue;
        }
        ledger.returned_contigs += kept.len();
        ledger.returned_bp += features.bin_size(&kept);
        for contig in &kept {
            pool.remove(contig);
        }
        bins.insert(bin_id, kept);
    }

    ledger.left_contigs = pool.len();
    ledger.left_bp = bases(features, &pool);
    *unbinned = sorted(pool);
    ledger
}
