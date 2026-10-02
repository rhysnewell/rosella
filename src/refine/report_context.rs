use std::collections::{BTreeMap, HashMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::embedding::knn::KnnGraph;
use crate::embedding::metrics::{AggregateMetric, Point};
use crate::quality::Scorer;
use crate::refine::audit::{neighbour_weight, share};
use crate::refine::bin_stats::centroid;
use crate::refine::owners::owners;

struct Profile {
    centre: Point,
}

impl Profile {
    fn of(features: &ContigFeatures, contigs: &[usize]) -> Option<Self> {
        if contigs.len() < 2 {
            return None;
        }
        Some(Self {
            centre: centroid(features, contigs),
        })
    }

    fn to(&self, metric: &AggregateMetric, point: &Point) -> f64 {
        metric.distance(point, &self.centre)
    }
}

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
    profiles: HashMap<usize, Profile>,
    families: HashMap<usize, HashSet<u32>>,
    carriers: HashMap<usize, HashMap<u32, usize>>,
}

impl Context {
    pub fn of(
        bins: &BTreeMap<usize, Vec<usize>>,
        features: &ContigFeatures<'_>,
        quality: &dyn Scorer,
    ) -> Self {
        let metric = AggregateMetric::new(features.n_samples() * 2, features.distance_settings());
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

        let profiles = members
            .iter()
            .filter_map(|(label, held)| {
                Profile::of(features, held).map(|profile| (*label, profile))
            })
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
            profiles,
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
            self.profiles
                .get(bin)
                .map_or(1.0, |profile| profile.to(&self.metric, &point))
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

    pub fn owned(&self) -> Vec<(usize, usize)> {
        let mut owned = self
            .owner
            .iter()
            .map(|(contig, label)| (*contig, *label))
            .collect::<Vec<_>>();
        owned.sort_unstable();
        owned
    }
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
