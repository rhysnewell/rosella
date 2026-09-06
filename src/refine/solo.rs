use std::collections::BTreeMap;

use crate::embedding::{features::ContigFeatures, metrics::AggregateMetric};
use crate::refine::bin_stats::{Centroid, centroid};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SoloPool {
    Alone,
    Majority,
    Long,
}

pub const SOLO_POOL_NAMES: [&str; 3] = ["alone", "majority", "long"];

impl SoloPool {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "alone" => Some(Self::Alone),
            "majority" => Some(Self::Majority),
            "long" => Some(Self::Long),
            _ => None,
        }
    }
}

/// Two strains at one depth look alike in every feature, and so do two halves of one genome.
/// Only scale separates them, and the run's own closed genomes say what genome-sized is here.
/// Mapped depth lets the partition fold closed genomes into fused bins, so `Alone` finds
/// nothing there; a contig holding most of its bin's bases is still that bin's genome.
pub fn floor(
    features: &ContigFeatures,
    bins: &BTreeMap<usize, Vec<usize>>,
    unbinned: &[usize],
    min_bin_size: usize,
    pool: SoloPool,
) -> Option<usize> {
    let mut alone = bins
        .values()
        .filter_map(|contigs| match pool {
            SoloPool::Alone => (contigs.len() == 1).then_some(vec![features.length(contigs[0])]),
            SoloPool::Majority => {
                let total = features.bin_size(contigs);
                let holders = contigs
                    .iter()
                    .map(|contig| features.length(*contig))
                    .filter(|length| contigs.len() == 1 || *length * 2 > total)
                    .collect::<Vec<_>>();
                (!holders.is_empty()).then_some(holders)
            }
            SoloPool::Long => Some(
                contigs
                    .iter()
                    .map(|contig| features.length(*contig))
                    .collect(),
            ),
        })
        .flatten()
        .chain(unbinned.iter().map(|contig| features.length(*contig)))
        .filter(|length| *length >= min_bin_size)
        .collect::<Vec<_>>();
    if alone.is_empty() {
        return None;
    }
    alone.sort_unstable();
    Some(alone[alone.len() / 2] / 2)
}

/// Sweeping every short contig into one leftover bin strands a genome's own fragments away
/// from the long contig just peeled out of their bin. Each is offered the piece it sits
/// nearest, with the leftover competing for it on the same terms.
fn scatter(
    features: &ContigFeatures,
    mut kept: Vec<Vec<usize>>,
    rest: Vec<usize>,
) -> Vec<Vec<usize>> {
    let anchors = kept
        .iter()
        .map(|piece| centroid(features, piece))
        .collect::<Vec<_>>();
    let metric = AggregateMetric::new(features.n_samples() * 2, features.distance_settings());
    let rows = features.rows(&rest);
    let mut sum = vec![0.0; rows[0].len()];
    let mut floors = 0.0;
    let mut total = 0.0;
    for (contig, row) in rest.iter().zip(&rows) {
        let weight = features.length(*contig) as f64;
        for (slot, value) in sum.iter_mut().zip(row) {
            *slot += value * weight;
        }
        floors += features.variance_floor(*contig) * weight;
        total += weight;
    }

    let mut leftover = Vec::new();
    for (contig, row) in rest.iter().zip(&rows) {
        let weight = features.length(*contig) as f64;
        let held = total - weight;
        let stay = (held > 0.0).then(|| {
            let others = Centroid {
                row: sum
                    .iter()
                    .zip(row)
                    .map(|(slot, value)| (slot - value * weight) / held)
                    .collect(),
                floor: (floors - features.variance_floor(*contig) * weight) / held,
            };
            metric.distance(
                row,
                &others.row,
                features.variance_floor(*contig),
                others.floor,
            )
        });
        let nearest = anchors
            .iter()
            .map(|anchor| {
                metric.distance(
                    row,
                    &anchor.row,
                    features.variance_floor(*contig),
                    anchor.floor,
                )
            })
            .enumerate()
            .min_by(|one, other| one.1.total_cmp(&other.1));
        match nearest {
            Some((piece, distance)) if stay.is_none_or(|stay| distance < stay) => {
                kept[piece].push(*contig)
            }
            _ => leftover.push(*contig),
        }
    }
    for piece in kept.iter_mut() {
        piece.sort_unstable();
    }
    if !leftover.is_empty() {
        kept.push(leftover);
    }
    kept
}

pub fn candidate(
    features: &ContigFeatures,
    indices: &[usize],
    floor: usize,
    scattered: bool,
) -> Option<Vec<Vec<usize>>> {
    let (alone, rest): (Vec<usize>, Vec<usize>) = indices
        .iter()
        .partition(|contig| features.length(**contig) >= floor);
    if alone.len() < 2 {
        return None;
    }
    let mut kept = match features.homology() {
        Some(homology) => homology.groups(&alone),
        None => alone.iter().map(|contig| vec![*contig]).collect(),
    };
    if kept.len() < 2 {
        return None;
    }
    if !rest.is_empty() {
        kept = if scattered {
            scatter(features, kept, rest)
        } else {
            kept.push(rest);
            kept
        };
    }
    Some(kept)
}
