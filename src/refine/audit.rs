use std::collections::{BTreeMap, HashMap, HashSet};

use rayon::prelude::*;

use crate::embedding::knn::KnnGraph;
use crate::refine::owners::owners;

/// A contig keeps its bin while the sequence around it agrees. The partition gives every contig
/// a label whatever its own evidence is worth, and nothing downstream asks again.
pub fn audit(
    bins: &mut BTreeMap<usize, Vec<usize>>,
    unbinned: &mut Vec<usize>,
    knn: &KnnGraph,
    lengths: &[usize],
) -> usize {
    let owner = owners(bins.iter().map(|(label, contigs)| (*label, contigs)));

    let mut evicted = owner
        .par_iter()
        .filter(|(contig, label)| {
            lengths[**contig] < crate::tuning::AUDIT_LENGTH
                && evict(**contig, **label, &owner, knn, lengths)
        })
        .map(|(contig, _)| *contig)
        .collect::<Vec<_>>();
    evicted.sort_unstable();

    let gone = evicted.iter().copied().collect::<HashSet<_>>();
    for members in bins.values_mut() {
        members.retain(|contig| !gone.contains(contig));
    }
    bins.retain(|_, members| !members.is_empty());
    let dropped = evicted.len();
    unbinned.append(&mut evicted);
    dropped
}

fn evict(
    contig: usize,
    label: usize,
    owner: &HashMap<usize, usize>,
    knn: &KnnGraph,
    lengths: &[usize],
) -> bool {
    let mut weights = neighbour_weight(contig, owner, knn, lengths);
    share(&mut weights, label).is_some_and(|share| share < crate::tuning::AUDIT_BAR)
}

/// Summed in bin order rather than the order the neighbours arrived, because a hasher seeded
/// per process would otherwise move the total by an ulp and flip a contig sitting on the bar.
pub fn share(weights: &mut [(usize, f64)], label: usize) -> Option<f64> {
    weights.sort_unstable_by_key(|(bin, _)| *bin);
    let total = weights.iter().map(|(_, weight)| *weight).sum::<f64>();
    if total <= 0.0 {
        return None;
    }
    let own = weights
        .iter()
        .find(|(bin, _)| *bin == label)
        .map_or(0.0, |(_, weight)| *weight);
    Some(own / total)
}

/// Weighted by neighbour length, so a contig's bin is judged by how much sequence backs it
/// rather than by how many neighbours it happens to have.
pub fn neighbour_weight(
    contig: usize,
    owner: &HashMap<usize, usize>,
    knn: &KnnGraph,
    lengths: &[usize],
) -> Vec<(usize, f64)> {
    let mut weights = Vec::<(usize, f64)>::new();
    if contig >= knn.n_points() {
        return weights;
    }
    for (neighbour, distance) in knn
        .indices
        .row(contig)
        .iter()
        .zip(knn.dists.row(contig).iter())
    {
        let neighbour = *neighbour as usize;
        if neighbour == contig {
            continue;
        }
        let Some(label) = owner.get(&neighbour) else {
            continue;
        };
        let weight = (1.0 - f64::from(*distance)).max(0.0) * lengths[neighbour] as f64;
        match weights.iter_mut().find(|(bin, _)| bin == label) {
            Some((_, held)) => *held += weight,
            None => weights.push((*label, weight)),
        }
    }
    weights
}
