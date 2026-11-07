use std::collections::HashMap;

use crate::embedding::knn::KnnGraph;

// Summed in bin order, so a contig's share does not move with the order its neighbours arrived in.
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

// Weighted by neighbour length, so a contig's bin is judged by how much sequence backs it
// rather than by how many neighbours it happens to have.
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
        crate::refine::owners::credit(&mut weights, *label, weight);
    }
    weights
}
