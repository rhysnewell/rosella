use std::collections::{BTreeMap, HashSet};

use rayon::prelude::*;

use crate::embedding::features::ContigFeatures;
use crate::refine::bar::{MIN_SPLIT_CONTIGS, levels};
use crate::refine::bin_stats::{AGGREGATE, BinStats, LevelSource, Thresholds, bin_stats};

/// Recruit, the outlier rescue and merge all add contigs to bins, and nothing takes one back
/// out short of a whole successful split. This is the other direction.
pub fn eject(
    features: &ContigFeatures,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    source: LevelSource,
    quantile: f64,
    factor: f64,
    min_bin_size: usize,
    seed: u64,
) -> Vec<usize> {
    let stats = bin_statistics(features, bins, seed);
    if stats.is_empty() {
        return Vec::new();
    }

    let thresholds = Thresholds::from_bins(
        stats
            .iter()
            .map(|(bin_id, stats)| (features.bin_size(&bins[bin_id]), stats)),
        source,
        quantile,
    );
    let bar = levels(&thresholds)[AGGREGATE] * factor;

    let mut ejected = Vec::new();
    for (bin_id, stats) in stats.iter() {
        let contigs = &bins[bin_id];
        let mut leaving = over_bar(contigs, stats, bar);
        if leaving.is_empty() {
            continue;
        }

        leaving.sort_by(|left, right| {
            stats.per_contig[right.1][AGGREGATE].total_cmp(&stats.per_contig[left.1][AGGREGATE])
        });

        let mut retained = features.bin_size(contigs);
        let mut taken = HashSet::new();
        for (contig, _) in leaving {
            let shrunk = retained - features.length(contig);
            if shrunk < min_bin_size {
                break;
            }
            retained = shrunk;
            taken.insert(contig);
        }

        if taken.is_empty() {
            continue;
        }
        bins.get_mut(bin_id)
            .expect("bin was read from this map")
            .retain(|contig| !taken.contains(contig));
        ejected.extend(taken);
    }

    ejected.sort_unstable();
    ejected
}

/// Aggregate alone, because the three column test in `misplaced_length` only has to nominate a
/// bin for re-clustering and fires on nearly every one. It is also the column recruit and
/// merge set their own bars on.
fn over_bar(contigs: &[usize], stats: &BinStats, bar: f64) -> Vec<(usize, usize)> {
    contigs
        .iter()
        .enumerate()
        .filter(|(position, _)| stats.per_contig[*position][AGGREGATE] > bar)
        .map(|(position, contig)| (*contig, position))
        .collect()
}

fn bin_statistics(
    features: &ContigFeatures,
    bins: &BTreeMap<usize, Vec<usize>>,
    seed: u64,
) -> BTreeMap<usize, BinStats> {
    bins.par_iter()
        .filter(|(_, contigs)| contigs.len() >= MIN_SPLIT_CONTIGS)
        .filter_map(|(bin_id, contigs)| {
            bin_stats(features, contigs, seed).map(|stats| (*bin_id, stats))
        })
        .collect()
}
