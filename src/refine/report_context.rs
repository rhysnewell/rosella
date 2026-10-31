use std::collections::{BTreeMap, HashMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::embedding::knn::KnnGraph;
use crate::embedding::metrics::{AggregateMetric, Point};
use crate::quality::Scorer;
use crate::refine::bin_stats::centroid;
use crate::refine::owners::owners;

// A chimera's centroid sits between the genomes it holds, so every member is equally far and
// its spread collapses. Standardising by that spread would have it claim every contig.
fn claim(taking: f64, taking_odds: f64, leaving: f64, leaving_odds: f64) -> f64 {
    let held = (1.0 - taking).clamp(0.0, 1.0) * taking_odds;
    let lost = (1.0 - leaving).clamp(0.0, 1.0) * leaving_odds;
    let total = held + lost;
    if total <= 0.0 { 0.0 } else { held / total }
}

fn wanted(families: &HashSet<u32>, held: &HashSet<u32>) -> f64 {
    let novel = families.difference(held).count();
    let duplicate = families.intersection(held).count();
    (1 + novel) as f64 / (1 + duplicate) as f64
}

// A family the bin holds only on this contig leaves with it.
fn needed(carried: &HashSet<u32>, carriers: &HashMap<u32, usize>, members: usize) -> f64 {
    if members <= 1 {
        return 1.0;
    }
    let sole = carried
        .iter()
        .filter(|family| carriers.get(family) == Some(&1))
        .count();
    (1 + sole) as f64
}

fn carriers(quality: &dyn Scorer, members: &[usize]) -> HashMap<u32, usize> {
    let mut carriers = HashMap::new();
    for contig in members {
        for family in quality.features(&[*contig]) {
            *carriers.entry(family).or_default() += 1;
        }
    }
    carriers
}

pub struct Context {
    pub metric: AggregateMetric,
    pub members: BTreeMap<usize, Vec<usize>>,
    owner: HashMap<usize, usize>,
    centres: HashMap<usize, Point>,
    families: HashMap<usize, HashSet<u32>>,
    carriers: HashMap<usize, HashMap<u32, usize>>,
}

impl Context {
    pub fn of(
        bins: &BTreeMap<usize, Vec<usize>>,
        features: &ContigFeatures<'_>,
        quality: &dyn Scorer,
    ) -> Self {
        let metric = AggregateMetric::new(
            features.n_samples() * 2,
            features.distance_settings().clone(),
        );
        let members = bins
            .iter()
            .map(|(label, contigs)| {
                (
                    *label,
                    crate::refine::ranking::sorted(contigs.iter().copied()),
                )
            })
            .collect::<BTreeMap<_, _>>();
        let owner = owners(members.iter().map(|(label, held)| (*label, held)));

        let centres = members
            .iter()
            .filter(|(_, held)| held.len() >= 2)
            .map(|(label, held)| (*label, centroid(features, held)))
            .collect();
        let families = members
            .iter()
            .map(|(label, held)| (*label, quality.features(held)))
            .collect();
        let carriers = members
            .iter()
            .map(|(label, held)| (*label, carriers(quality, held)))
            .collect();

        Self {
            metric,
            members,
            owner,
            centres,
            families,
            carriers,
        }
    }

    pub fn labels(&self) -> Vec<usize> {
        self.members.keys().copied().collect()
    }

    pub fn rivals(
        &self,
        contig: usize,
        label: usize,
        inputs: &Inputs<'_>,
    ) -> Option<(f64, Vec<Rival>)> {
        let mut weights = neighbour_weight(contig, &self.owner, inputs.knn, inputs.lengths);
        let own = share(&mut weights, label)?;
        let total = weights.iter().map(|(_, weight)| *weight).sum::<f64>();
        let point = inputs.features.point(contig);
        let carried = inputs.quality.features(&[contig]);
        let distance = |bin: &usize| {
            self.centres
                .get(bin)
                .map_or(1.0, |centre| self.metric.distance(&point, centre))
        };
        let leaving = distance(&label);
        let leaving_odds = match (self.members.get(&label), self.carriers.get(&label)) {
            (Some(donor), Some(carriers)) => needed(&carried, carriers, donor.len()),
            _ => 1.0,
        };
        let rivals = weights
            .iter()
            .filter(|(bin, _)| *bin != label)
            .map(|(bin, weight)| {
                let odds = self
                    .families
                    .get(bin)
                    .map_or(1.0, |families| wanted(&carried, families));
                Rival {
                    bin: *bin,
                    share: weight / total,
                    claim: claim(distance(bin), odds, leaving, leaving_odds),
                }
            })
            .collect();
        Some((own, rivals))
    }
}

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

pub struct Inputs<'a> {
    pub features: &'a ContigFeatures<'a>,
    pub quality: &'a dyn Scorer,
    pub knn: &'a KnnGraph,
    pub lengths: &'a [usize],
    pub names: &'a [String],
}

pub struct Rival {
    pub bin: usize,
    pub share: f64,
    pub claim: f64,
}
