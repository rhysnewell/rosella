use std::collections::{BTreeMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::embedding::metrics::AggregateMetric;
use crate::quality::Scorer;
use crate::refine::bin_stats::centroid;

#[derive(Debug, Clone, Copy)]
pub struct RecruitSettings {
    pub completeness: f64,
    pub contamination: f64,
    pub margin: f64,
    pub min_bin_size: usize,
    pub worth: f64,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct RecruitLedger {
    pub near_bar: usize,
    pub donor_bins: usize,
    pub donor_contigs: usize,
    pub offered: usize,
    pub taken: usize,
    pub taken_bp: usize,
    pub lifted: usize,
}

impl std::fmt::Display for RecruitLedger {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "{} bins within the margin drew on {} contigs from {} bins under the floor; \
             offered {}, took {} holding {} bp, lifting {} bins over the bar",
            self.near_bar,
            self.donor_contigs,
            self.donor_bins,
            self.offered,
            self.taken,
            self.taken_bp,
            self.lifted
        )
    }
}

/// A bin under the size floor is discarded when the bins are written, so offering its contigs
/// to a bin a few points short of the bar risks nothing and is the only place those bases can
/// still land.
pub fn recruit(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    settings: RecruitSettings,
) -> RecruitLedger {
    let mut ledger = RecruitLedger::default();
    let floor = settings.completeness - settings.margin;

    let mut donors: Vec<usize> = Vec::new();
    let mut near: Vec<(f64, usize)> = Vec::new();
    for (id, contigs) in bins.iter() {
        if features.bin_size(contigs) < settings.min_bin_size {
            ledger.donor_bins += 1;
            donors.extend(contigs.iter().copied());
            continue;
        }
        let held = quality.score(contigs);
        if held.contamination > settings.contamination || held.completeness >= settings.completeness
        {
            continue;
        }
        if held.completeness >= floor {
            near.push((held.completeness, *id));
        }
    }
    ledger.near_bar = near.len();
    ledger.donor_contigs = donors.len();
    if near.is_empty() || donors.is_empty() {
        return ledger;
    }

    donors.sort_unstable();
    near.sort_by(|one, other| other.0.total_cmp(&one.0));

    let metric = AggregateMetric::new(features.n_samples() * 2, features.distance_settings());
    let mut claimed: HashSet<usize> = HashSet::new();
    let mut moved: Vec<usize> = Vec::new();

    for (_, id) in near {
        let mut contigs = bins[&id].clone();
        let mut held = quality.score(&contigs);
        let mut families = quality.features(&contigs);
        let centre = centroid(features, &contigs);

        let mut offers = donors
            .iter()
            .filter(|contig| !claimed.contains(contig))
            .filter(|contig| !quality.features(&[**contig]).is_subset(&families))
            .map(|contig| (metric.distance(&row(features, *contig), &centre.row, crate::embedding::metrics::MIN_VAR, centre.floor), *contig))
            .collect::<Vec<_>>();
        offers.sort_by(|one, other| one.0.total_cmp(&other.0));

        let mut took = false;
        for (_, contig) in offers {
            ledger.offered += 1;
            let mut trial = contigs.clone();
            trial.push(contig);
            trial.sort_unstable();
            let scored = quality.score(&trial);
            if scored.contamination > settings.contamination
                || scored.score(settings.worth) <= held.score(settings.worth)
            {
                continue;
            }
            contigs = trial;
            held = scored;
            families = quality.features(&contigs);
            claimed.insert(contig);
            moved.push(contig);
            ledger.taken += 1;
            ledger.taken_bp += features.length(contig);
            took = true;
            if held.completeness >= settings.completeness {
                break;
            }
        }
        if took {
            if held.completeness >= settings.completeness {
                ledger.lifted += 1;
            }
            bins.insert(id, contigs);
        }
    }

    if !moved.is_empty() {
        let taken = moved.into_iter().collect::<HashSet<_>>();
        for contigs in bins.values_mut() {
            if features.bin_size(contigs) < settings.min_bin_size {
                contigs.retain(|contig| !taken.contains(contig));
            }
        }
        bins.retain(|_, contigs| !contigs.is_empty());
    }
    ledger
}

fn row(features: &ContigFeatures, contig: usize) -> Vec<f64> {
    let mut held = Vec::with_capacity(features.n_samples() * 2 + features.tnf_row(contig).len());
    held.extend_from_slice(features.coverage_row(contig));
    held.extend_from_slice(features.tnf_row(contig));
    held
}
