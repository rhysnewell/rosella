use std::collections::{BTreeMap, HashMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::embedding::metrics::AggregateMetric;
use crate::quality::Scorer;
use crate::refine::bin_stats::centroid;

pub struct Profile {
    centre: Vec<f64>,
    floor: f64,
}

impl Profile {
    pub fn of(features: &ContigFeatures, contigs: &[usize]) -> Option<Self> {
        if contigs.len() < 2 {
            return None;
        }
        let centre = centroid(features, contigs);
        Some(Self {
            centre: centre.row,
            floor: centre.floor,
        })
    }

    pub fn to(&self, metric: &AggregateMetric, row: &[f64], floor: f64) -> f64 {
        metric.distance(row, &self.centre, floor, self.floor)
    }
}

// A chimera's centroid sits between the genomes it holds, so every member is equally far and
// its spread collapses. Standardising by that spread would have it claim every contig.
pub fn claim(taking: f64, taking_odds: f64, leaving: f64, leaving_odds: f64) -> f64 {
    let held = (1.0 - taking).clamp(0.0, 1.0) * taking_odds;
    let lost = (1.0 - leaving).clamp(0.0, 1.0) * leaving_odds;
    let total = held + lost;
    if total <= 0.0 { 0.0 } else { held / total }
}

pub fn wanted(families: &HashSet<u32>, held: &HashSet<u32>) -> f64 {
    let novel = families.difference(held).count();
    let duplicate = families.intersection(held).count();
    (1 + novel) as f64 / (1 + duplicate) as f64
}

pub fn needed(quality: &dyn Scorer, donor: &[usize], contig: usize) -> f64 {
    let left = donor
        .iter()
        .copied()
        .filter(|other| *other != contig)
        .collect::<Vec<_>>();
    if left.is_empty() {
        return 1.0;
    }
    let sole = quality
        .features(donor)
        .difference(&quality.features(&left))
        .count();
    (1 + sole) as f64
}

pub struct Context {
    pub metric: AggregateMetric,
    pub members: HashMap<usize, Vec<usize>>,
    pub owner: HashMap<usize, usize>,
    pub profiles: HashMap<usize, Profile>,
    pub families: HashMap<usize, HashSet<u32>>,
}

impl Context {
    pub fn of(
        bins: &BTreeMap<usize, Vec<usize>>,
        features: &ContigFeatures<'_>,
        quality: &dyn Scorer,
    ) -> Self {
        let metric = AggregateMetric::new(features.n_samples() * 2, features.distance_settings());
        let mut members: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut owner: HashMap<usize, usize> = HashMap::new();
        for (label, contigs) in bins {
            let mut held = contigs.clone();
            held.sort_unstable();
            for contig in &held {
                owner.insert(*contig, *label);
            }
            members.insert(*label, held);
        }

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

        Self {
            metric,
            members,
            owner,
            profiles,
            families,
        }
    }

    pub fn labels(&self) -> Vec<usize> {
        let mut labels = self.members.keys().copied().collect::<Vec<_>>();
        labels.sort_unstable();
        labels
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
