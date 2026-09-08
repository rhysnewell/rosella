use std::collections::{BTreeMap, HashMap, HashSet};

use rayon::prelude::*;

use crate::embedding::features::ContigFeatures;

pub const DEFAULT_BAR: f64 = 0.05;
pub const DEFAULT_LINK: f64 = 0.5;
pub const DEFAULT_MIN_HASHES: usize = 20;

/// A k-mer this widely held is a low complexity repeat, and pairing it costs more than it says.
pub const MAX_SPREAD: usize = 8;

#[derive(Debug, Clone, Copy)]
pub struct DuplicationSettings {
    pub bar: f64,
    pub link: f64,
    pub min_hashes: usize,
}

impl Default for DuplicationSettings {
    fn default() -> Self {
        Self {
            bar: DEFAULT_BAR,
            link: DEFAULT_LINK,
            min_hashes: DEFAULT_MIN_HASHES,
        }
    }
}

/// Containment is symmetric evidence, so it says two contigs share sequence and not which one
/// is the intruder. Mass breaks the tie: the side the rest of the bin swallows is the smaller.
pub fn intruders(
    features: &ContigFeatures,
    contigs: &[usize],
    settings: DuplicationSettings,
) -> Vec<(usize, f64)> {
    let Some(sketches) = features.sketches() else {
        return Vec::new();
    };
    if contigs.len() < 2 || sketches.duplication(contigs).unwrap_or(0.0) < settings.bar {
        return Vec::new();
    }

    let mut owners: HashMap<u64, Vec<u32>> = HashMap::new();
    for (position, contig) in contigs.iter().enumerate() {
        for hash in sketches.hashes(*contig) {
            owners.entry(*hash).or_default().push(position as u32);
        }
    }

    let mut reached = vec![0u32; contigs.len()];
    let mut shared: HashMap<(u32, u32), u32> = HashMap::new();
    for holders in owners.values() {
        if holders.len() < 2 || holders.len() > MAX_SPREAD {
            continue;
        }
        for position in holders {
            reached[*position as usize] += 1;
        }
        for left in 0..holders.len() {
            for right in left + 1..holders.len() {
                *shared.entry((holders[left], holders[right])).or_default() += 1;
                *shared.entry((holders[right], holders[left])).or_default() += 1;
            }
        }
    }

    let mut leaving = Vec::new();
    for (position, contig) in contigs.iter().enumerate() {
        let held = sketches.hashes(*contig).len();
        if held < settings.min_hashes {
            continue;
        }
        let containment = f64::from(reached[position]) / held as f64;
        if containment < settings.link {
            continue;
        }
        let partner_bases: usize = contigs
            .iter()
            .enumerate()
            .filter(|(other, _)| *other != position)
            .filter(|(other, _)| {
                let count = shared
                    .get(&(position as u32, *other as u32))
                    .copied()
                    .unwrap_or(0);
                f64::from(count) / held as f64 >= settings.link
            })
            .map(|(_, partner)| features.length(*partner))
            .sum();
        if partner_bases > features.length(*contig) {
            leaving.push((*contig, containment));
        }
    }
    leaving
}

/// The distance arm asks whether a contig sits apart from its binmates. This one asks whether
/// the bin holds the same sequence twice, which is what two fused strains look like.
pub fn eject_duplicated(
    features: &ContigFeatures,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    settings: DuplicationSettings,
    min_bin_size: usize,
) -> Vec<usize> {
    if features.sketches().is_none() {
        return Vec::new();
    }

    let proposals = bins
        .par_iter()
        .map(|(bin_id, contigs)| (*bin_id, intruders(features, contigs, settings)))
        .filter(|(_, leaving)| !leaving.is_empty())
        .collect::<Vec<_>>();

    let mut ejected = Vec::new();
    for (bin_id, mut leaving) in proposals {
        leaving.sort_by(|left, right| right.1.total_cmp(&left.1));
        let contigs = &bins[&bin_id];
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
        bins.get_mut(&bin_id)
            .expect("bin was read from this map")
            .retain(|contig| !taken.contains(contig));
        ejected.extend(taken);
    }

    ejected.sort_unstable();
    ejected
}
