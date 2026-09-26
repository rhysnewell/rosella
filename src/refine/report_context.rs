use std::collections::{BTreeMap, HashMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::embedding::metrics::AggregateMetric;
use crate::quality::Scorer;
use crate::refine::recruit::Profile;

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
                Profile::of(features, &metric, held).map(|profile| (*label, profile))
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
