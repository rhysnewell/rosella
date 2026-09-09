use std::collections::{BTreeMap, HashMap, HashSet};

use anyhow::Result;
use log::{info, warn};

use crate::clustering::clusterer::{Partitioning, placed_once};
use crate::embedding::features::ContigFeatures;
use crate::embedding::knn::KnnGraph;
use crate::quality::ContigQuality;
use crate::refine::rung::{Bars, Rung, Verdict, judge, over_bar};
use crate::refine::select::ranked;

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

#[derive(Debug, Clone, Copy)]
pub struct DissolveSettings {
    pub bars: Bars,
    pub genome_floor: Option<usize>,
    pub min_contigs: usize,
    pub rounds: usize,
    pub passes: usize,
    pub n_neighbours: usize,
    pub reuse: bool,
}

/// What the pool took, what it refused and where the refusals went, in contigs and bases.
/// The bin count alone cannot say whether an arm has headroom left.
#[derive(Debug, Default, Clone, Copy)]
pub struct DissolveLedger {
    pub dissolved_small: usize,
    pub dissolved_duplicated: usize,
    pub dissolved_clean: usize,
    pub dissolved_bp: usize,
    pub pool_contigs: usize,
    pub pool_bp: usize,
    pub rounds: usize,
    pub passes: usize,
    pub rung: usize,
    pub proposed: usize,
    pub proposed_composition: usize,
    pub refused_small: usize,
    pub refused_incomplete: usize,
    pub refused_contaminated: usize,
    pub refused_duplicated: usize,
    pub refused_worse: usize,
    pub noise: usize,
    pub promoted: usize,
    pub adopted_contigs: usize,
    pub adopted_bp: usize,
    pub returned_contigs: usize,
    pub returned_bp: usize,
    pub emptied: usize,
    pub left_contigs: usize,
    pub left_bp: usize,
}

impl std::fmt::Display for DissolveLedger {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "dissolved {} small, {} duplicated and {} clean bins holding {} bp; pool {} contigs \
             {} bp; {} rounds over {} passes ending at rung {}; proposed {} clusters, refused {} small, {} \
             incomplete, {} contaminated, {} duplicated and {} no better, {} noise; {} of the \
             proposals came only from composition; promoted {} \
             bins adopting {} contigs {} bp; returned {} contigs {} bp, emptied {} bins; left {} \
             contigs {} bp unbinned",
            self.dissolved_small,
            self.dissolved_duplicated,
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
            self.refused_duplicated,
            self.refused_worse,
            self.noise,
            self.proposed_composition,
            self.promoted,
            self.adopted_contigs,
            self.adopted_bp,
            self.returned_contigs,
            self.returned_bp,
            self.emptied,
            self.left_contigs,
            self.left_bp
        )
    }
}

/// A piece smaller than the run's own genome scale is a shard of one, and relaxing that to let
/// more of the pool through cost more bins than it recovered.
fn floor_for(settings: DissolveSettings) -> usize {
    settings
        .genome_floor
        .unwrap_or(settings.bars.min_bin_size)
        .max(settings.bars.min_bin_size)
}

fn sorted(contigs: HashSet<usize>) -> Vec<usize> {
    let mut contigs = contigs.into_iter().collect::<Vec<_>>();
    contigs.sort_unstable();
    contigs
}

fn bases(features: &ContigFeatures, contigs: &HashSet<usize>) -> usize {
    contigs.iter().map(|contig| features.length(*contig)).sum()
}

/// Every bin goes back in, because a bin that survived the earlier stages is still only what
/// composition and coverage could group, and the gene families judge it on a different axis.
fn dissolving(
    features: &ContigFeatures,
    bins: &BTreeMap<usize, Vec<usize>>,
    top: usize,
    settings: DissolveSettings,
    ledger: &mut DissolveLedger,
) -> Vec<(usize, Vec<usize>)> {
    let mut dissolving = Vec::new();
    for (bin_id, contigs) in bins.iter() {
        let small = features.bin_size(contigs) < top;
        let duplicated = !small && over_bar(features, contigs, settings.bars.duplication_bar);
        if small {
            ledger.dissolved_small += 1;
        } else if duplicated {
            ledger.dissolved_duplicated += 1;
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
    quality: Option<&'a ContigQuality>,
    origin: HashMap<usize, usize>,
    held: HashMap<usize, (f64, usize)>,
}

impl Pot<'_> {
    pub fn scored(&self) -> bool {
        self.quality.is_some()
    }

    pub fn worth(&self, contigs: &[usize]) -> f64 {
        match self.quality {
            Some(quality) => quality.score(contigs).score(),
            None => self.features.bin_size(contigs) as f64,
        }
    }

    pub fn judge(&self, contigs: &[usize], rung: Rung) -> Verdict {
        judge(self.features, self.quality, contigs, rung)
    }

    /// A cluster that takes the greater part of a bin has to be the better bin, or the loop
    /// trades a whole genome for a piece of one.
    pub fn improves(&self, contigs: &[usize]) -> bool {
        let mut taken: HashMap<usize, usize> = HashMap::new();
        for contig in contigs {
            if let Some(bin) = self.origin.get(contig) {
                *taken.entry(*bin).or_default() += self.features.length(*contig);
            }
        }
        let Some(quality) = self.quality else {
            return true;
        };
        let candidate = quality.score(contigs).score();
        taken
            .into_iter()
            .all(|(bin, bases)| match self.held.get(&bin) {
                Some((score, whole)) if bases * 2 >= *whole => candidate > *score,
                _ => true,
            })
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
    info!(
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
pub fn dissolve(
    features: &ContigFeatures,
    quality: Option<&ContigQuality>,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    unbinned: &mut Vec<usize>,
    settings: DissolveSettings,
    oracle: &[Vec<usize>],
    neighbours: impl Fn(&HashSet<usize>, usize, PoolView) -> Result<(KnnGraph, Vec<usize>)>,
    partition: impl Fn(&KnnGraph, &[usize], RoundParams) -> Result<Vec<Partitioning>>,
) -> DissolveLedger {
    let mut ledger = DissolveLedger::default();
    let top = floor_for(settings);
    let dissolved = dissolving(features, bins, top, settings, &mut ledger);

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
    let pot = Pot {
        features,
        quality,
        origin: dissolved
            .iter()
            .flat_map(|(bin_id, contigs)| contigs.iter().map(|contig| (*contig, *bin_id)))
            .collect(),
        held: match quality {
            Some(quality) => dissolved
                .iter()
                .map(|(bin_id, contigs)| {
                    (
                        *bin_id,
                        (quality.score(contigs).score(), features.bin_size(contigs)),
                    )
                })
                .collect(),
            None => HashMap::new(),
        },
    };
    let mut promoted = ranked(
        &pot,
        &mut pool,
        settings,
        oracle,
        top,
        &mut ledger,
        neighbours,
        partition,
    );
    if !oracle.is_empty() {
        report_oracle(features, &handed, oracle, &promoted);
    }
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

    // Whatever the pool did not claim goes back to the bin it came from, or a long contig the
    // writer would stand alone leaves as a singleton for no gain.
    for (bin_id, contigs) in dissolved {
        let kept = contigs
            .iter()
            .copied()
            .filter(|contig| !claimed.contains(contig))
            .collect::<Vec<_>>();
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
