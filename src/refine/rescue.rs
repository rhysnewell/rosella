use std::collections::{BTreeMap, HashSet};

use anyhow::Result;
use log::warn;

use crate::clustering::clusterer::HDBSCANResult;
use crate::embedding::features::ContigFeatures;

#[derive(Debug, Clone, Copy)]
pub struct RescueSettings {
    pub min_bin_size: usize,
    pub genome_floor: Option<usize>,
    pub duplication_bar: f64,
    pub min_contigs: usize,
}

/// A piece smaller than the run's own genome scale is a shard of one, and relaxing that to let
/// more of the pool through cost more bins than it recovered.
fn floor_for(settings: RescueSettings) -> usize {
    settings
        .genome_floor
        .unwrap_or(settings.min_bin_size)
        .max(settings.min_bin_size)
}

fn accepted(features: &ContigFeatures, contigs: &[usize], floor: usize, bar: f64) -> bool {
    features.bin_size(contigs) >= floor
        && features
            .sketches()
            .and_then(|sketches| sketches.duplication(contigs))
            .is_none_or(|duplication| duplication <= bar)
}

/// The eject arm has already stripped what it could before this runs, so a bin still over the
/// bar is one it could not fix, which is the bin least worth trusting as it stands.
fn unsure(features: &ContigFeatures, contigs: &[usize], bar: f64) -> bool {
    features
        .sketches()
        .and_then(|sketches| sketches.duplication(contigs))
        .is_some_and(|duplication| duplication > bar)
}

fn sorted(contigs: HashSet<usize>) -> Vec<usize> {
    let mut contigs = contigs.into_iter().collect::<Vec<_>>();
    contigs.sort_unstable();
    contigs
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
) -> usize {
    let top = floor_for(settings);
    let dissolving = bins
        .iter()
        .filter(|(_, contigs)| {
            features.bin_size(contigs) < top
                || unsure(features, contigs, settings.duplication_bar)
        })
        .map(|(bin_id, contigs)| (*bin_id, contigs.clone()))
        .collect::<Vec<_>>();

    let mut pool = unbinned.iter().copied().collect::<HashSet<_>>();
    for (_, contigs) in &dissolving {
        pool.extend(contigs.iter().copied());
    }

    if pool.len() < settings.min_contigs {
        return 0;
    }
    let result = match partition(&pool) {
        Ok(result) => result,
        Err(error) => {
            warn!("Could not re-embed the rescue pool: {error}");
            return 0;
        }
    };
    let promoted = result
        .cluster_map
        .into_values()
        .map(sorted)
        .filter(|contigs| accepted(features, contigs, top, settings.duplication_bar))
        .collect::<Vec<_>>();

    if promoted.is_empty() {
        return 0;
    }

    let next_bin_id = bins.keys().max().map_or(0, |id| id + 1);
    let mut claimed = HashSet::new();
    let taken = promoted.len();
    for (offset, contigs) in promoted.into_iter().enumerate() {
        for contig in &contigs {
            pool.remove(contig);
        }
        claimed.extend(contigs.iter().copied());
        bins.insert(next_bin_id + offset, contigs);
    }

    // Whatever the pool did not claim goes back to the bin it came from, or a long contig the
    // writer would stand alone leaves as a singleton for no gain.
    for (bin_id, contigs) in dissolving {
        let kept = contigs
            .iter()
            .copied()
            .filter(|contig| !claimed.contains(contig))
            .collect::<Vec<_>>();
        if kept.is_empty() {
            bins.remove(&bin_id);
        } else {
            for contig in &kept {
                pool.remove(contig);
            }
            bins.insert(bin_id, kept);
        }
    }
    *unbinned = sorted(pool);
    taken
}
