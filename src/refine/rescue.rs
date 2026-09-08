use std::collections::{BTreeMap, HashSet};

use anyhow::Result;
use log::warn;

use crate::clustering::clusterer::{HDBSCANResult, placed_once};
use crate::embedding::features::ContigFeatures;

#[derive(Debug, Clone, Copy)]
pub struct RescueSettings {
    pub min_bin_size: usize,
    pub genome_floor: Option<usize>,
    pub duplication_bar: f64,
    pub min_contigs: usize,
}

/// What the pool took, what it refused and where the refusals went, in contigs and bases.
/// The bin count alone cannot say whether an arm has headroom left.
#[derive(Debug, Default, Clone, Copy)]
pub struct RescueLedger {
    pub dissolved_small: usize,
    pub dissolved_duplicated: usize,
    pub dissolved_bp: usize,
    pub pool_contigs: usize,
    pub pool_bp: usize,
    pub proposed: usize,
    pub refused_small: usize,
    pub refused_duplicated: usize,
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

impl std::fmt::Display for RescueLedger {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "dissolved {} small and {} duplicated bins holding {} bp; pool {} contigs {} bp; \
             proposed {} clusters, refused {} small and {} duplicated, {} noise; \
             promoted {} bins adopting {} contigs {} bp; returned {} contigs {} bp, \
             emptied {} bins; left {} contigs {} bp unbinned",
            self.dissolved_small,
            self.dissolved_duplicated,
            self.dissolved_bp,
            self.pool_contigs,
            self.pool_bp,
            self.proposed,
            self.refused_small,
            self.refused_duplicated,
            self.noise,
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
fn floor_for(settings: RescueSettings) -> usize {
    settings
        .genome_floor
        .unwrap_or(settings.min_bin_size)
        .max(settings.min_bin_size)
}

fn over_bar(features: &ContigFeatures, contigs: &[usize], bar: f64) -> bool {
    features
        .sketches()
        .and_then(|sketches| sketches.duplication(contigs))
        .is_some_and(|duplication| duplication > bar)
}

enum Verdict {
    Adopt,
    TooSmall,
    Duplicated,
}

fn judge(features: &ContigFeatures, contigs: &[usize], floor: usize, bar: f64) -> Verdict {
    if features.bin_size(contigs) < floor {
        Verdict::TooSmall
    } else if over_bar(features, contigs, bar) {
        Verdict::Duplicated
    } else {
        Verdict::Adopt
    }
}

fn sorted(contigs: HashSet<usize>) -> Vec<usize> {
    let mut contigs = contigs.into_iter().collect::<Vec<_>>();
    contigs.sort_unstable();
    contigs
}

fn bases(features: &ContigFeatures, contigs: &HashSet<usize>) -> usize {
    contigs.iter().map(|contig| features.length(*contig)).sum()
}

/// The eject arm has already stripped what it could before this runs, so a bin still over the
/// bar is one it could not fix, which is the bin least worth trusting as it stands.
fn dissolving(
    features: &ContigFeatures,
    bins: &BTreeMap<usize, Vec<usize>>,
    top: usize,
    settings: RescueSettings,
    ledger: &mut RescueLedger,
) -> Vec<(usize, Vec<usize>)> {
    let mut dissolving = Vec::new();
    for (bin_id, contigs) in bins.iter() {
        let small = features.bin_size(contigs) < top;
        if !small && !over_bar(features, contigs, settings.duplication_bar) {
            continue;
        }
        if small {
            ledger.dissolved_small += 1;
        } else {
            ledger.dissolved_duplicated += 1;
        }
        ledger.dissolved_bp += features.bin_size(contigs);
        dissolving.push((*bin_id, contigs.clone()));
    }
    dissolving
}

/// Bins under the genome floor, and bins the sketch says hold their own sequence twice, go back
/// in the pot with the unbinned and are embedded again without the bins that already left, which
/// is the one thing re-cutting inside a bin cannot do.
pub fn rescue(
    features: &ContigFeatures,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    unbinned: &mut Vec<usize>,
    settings: RescueSettings,
    partition: impl Fn(&HashSet<usize>) -> Result<HDBSCANResult>,
) -> RescueLedger {
    let mut ledger = RescueLedger::default();
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
    let result = match partition(&pool) {
        Ok(result) => result,
        Err(error) => {
            warn!("Could not re-embed the rescue pool: {error}");
            return ledger;
        }
    };

    ledger.noise = result.outliers.len();
    let mut promoted = Vec::new();
    for contigs in result.cluster_map.into_values().map(sorted) {
        ledger.proposed += 1;
        match judge(features, &contigs, top, settings.duplication_bar) {
            Verdict::Adopt => promoted.push(contigs),
            Verdict::TooSmall => ledger.refused_small += 1,
            Verdict::Duplicated => ledger.refused_duplicated += 1,
        }
    }
    if promoted.is_empty() {
        return ledger;
    }
    // Hash order would give the same partition different bin names on every run.
    promoted.sort_unstable_by_key(|contigs| contigs.first().copied().unwrap_or(usize::MAX));

    if let Err(error) = placed_once(promoted.iter().flatten().copied(), &pool) {
        warn!("The rescue pool did not come back partitioned, leaving the bins alone: {error}");
        return ledger;
    }

    let next_bin_id = bins.keys().max().map_or(0, |id| id + 1);
    let mut claimed = HashSet::new();
    ledger.promoted = promoted.len();
    for (offset, contigs) in promoted.into_iter().enumerate() {
        ledger.adopted_contigs += contigs.len();
        ledger.adopted_bp += features.bin_size(&contigs);
        for contig in &contigs {
            pool.remove(contig);
        }
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
