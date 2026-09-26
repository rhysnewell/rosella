use std::collections::{BTreeMap, HashMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::embedding::knn::KnnGraph;
use crate::embedding::metrics::AggregateMetric;
use crate::quality::{Quality, Scorer};
use crate::refine::bin_stats::centroid;

pub const DEFAULT_FLOOR: f64 = 0.78;
pub const DEFAULT_CONFIDENCE: f64 = 0.65;

#[derive(Debug, Default, Clone, Copy)]
pub struct RecruitLedger {
    pub bins: usize,
    pub receivers: usize,
    pub offered: usize,
    pub taken: usize,
    pub taken_bp: usize,
    pub emptied: usize,
    pub passes: usize,
}

impl std::fmt::Display for RecruitLedger {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "{} bins over {} passes, {} recruiting; offered {} neighbours, took {} contigs \
             over {} bp, emptied {}",
            self.bins,
            self.passes,
            self.receivers,
            self.offered,
            self.taken,
            self.taken_bp,
            self.emptied
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RecruitSettings {
    pub floor: f64,
    pub confidence: f64,
    pub completeness: f64,
    pub contamination: f64,
    pub max_bin_size: usize,
    pub passes: usize,
}

pub struct Profile {
    centre: Vec<f64>,
    floor: f64,
    mean: f64,
    deviation: f64,
}

impl Profile {
    pub fn of(
        features: &ContigFeatures,
        metric: &AggregateMetric,
        contigs: &[usize],
    ) -> Option<Self> {
        if contigs.len() < 2 {
            return None;
        }
        let centre = centroid(features, contigs);
        let floors = features.floors(contigs);
        let spread = features
            .rows(contigs)
            .iter()
            .zip(&floors)
            .map(|(row, floor)| metric.distance(row, &centre.row, *floor, centre.floor))
            .collect::<Vec<_>>();
        let mean = spread.iter().sum::<f64>() / spread.len() as f64;
        let variance = spread
            .iter()
            .map(|distance| (distance - mean) * (distance - mean))
            .sum::<f64>()
            / spread.len() as f64;
        Some(Self {
            centre: centre.row,
            floor: centre.floor,
            mean,
            deviation: variance.sqrt(),
        })
    }

    pub fn to(&self, metric: &AggregateMetric, row: &[f64], floor: f64) -> f64 {
        metric.distance(row, &self.centre, floor, self.floor)
    }

    pub fn strays(&self, distance: f64) -> bool {
        distance > self.mean + crate::tuning::PEEL_SIGMA * self.deviation
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

fn owner(bins: &BTreeMap<usize, Vec<usize>>) -> HashMap<usize, usize> {
    let mut of = HashMap::new();
    for (id, contigs) in bins {
        for contig in contigs {
            of.insert(*contig, *id);
        }
    }
    of
}

fn neighbours(knn: &KnnGraph, contigs: &[usize]) -> Vec<usize> {
    let members = contigs.iter().copied().collect::<HashSet<_>>();
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for contig in contigs {
        if *contig >= knn.n_points() {
            continue;
        }
        for neighbour in knn.indices.row(*contig) {
            let neighbour = *neighbour as usize;
            if !members.contains(&neighbour) && seen.insert(neighbour) {
                out.push(neighbour);
            }
        }
    }
    out.sort_unstable();
    out
}

fn whole(held: Quality, settings: RecruitSettings) -> bool {
    held.clears(crate::quality::Bars {
        completeness: settings.completeness,
        contamination: settings.contamination,
    })
}

fn with(contigs: &[usize], contig: usize) -> Vec<usize> {
    let mut joined = Vec::with_capacity(contigs.len() + 1);
    joined.extend_from_slice(contigs);
    joined.push(contig);
    joined.sort_unstable();
    joined
}

fn pass(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    knn: &KnnGraph,
    metric: &AggregateMetric,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    settings: RecruitSettings,
    ledger: &mut RecruitLedger,
) -> usize {
    // Stale as the pass runs, because a move changes both bins: the passes refresh the
    // profiles rather than every move paying for a centroid.
    let profiles = bins
        .iter()
        .filter_map(|(id, contigs)| Profile::of(features, metric, contigs).map(|held| (*id, held)))
        .collect::<BTreeMap<_, _>>();
    let scores = bins
        .iter()
        .map(|(id, contigs)| (*id, quality.score(contigs)))
        .collect::<BTreeMap<_, _>>();

    let mut receivers = bins
        .keys()
        .copied()
        .filter(|id| {
            profiles.contains_key(id)
                && scores[id].completeness >= settings.floor
                && scores[id].contamination <= settings.contamination
        })
        .collect::<Vec<_>>();
    receivers.sort_by(|left, right| {
        scores[right]
            .completeness
            .total_cmp(&scores[left].completeness)
            .then(left.cmp(right))
    });
    ledger.receivers = ledger.receivers.max(receivers.len());

    let mut of = owner(bins);
    let mut moved = 0;
    for id in receivers {
        let mut contigs = bins[&id].clone();
        if features.bin_size(&contigs) >= settings.max_bin_size {
            continue;
        }
        let taking = &profiles[&id];

        let candidates = neighbours(knn, &contigs)
            .into_iter()
            .filter(|contig| of.get(contig).is_some_and(|from| *from != id))
            .collect::<Vec<_>>();
        let floors = features.floors(&candidates);
        let rows = features.rows(&candidates);

        let families = quality.features(&contigs);
        let mut shortlist = Vec::new();
        for ((contig, row), floor) in candidates.iter().zip(&rows).zip(&floors) {
            let distance = taking.to(metric, row, *floor);
            if taking.strays(distance) {
                continue;
            }
            let losing = of
                .get(contig)
                .and_then(|from| profiles.get(from))
                .map_or(1.0, |leaving| leaving.to(metric, row, *floor));
            let odds = wanted(&quality.features(&[*contig]), &families);
            if claim(distance, odds, losing, 1.0) >= settings.confidence {
                shortlist.push((distance, losing, odds, *contig));
            }
        }
        ledger.offered += shortlist.len();
        shortlist.sort_by(|left, right| {
            claim(right.0, right.2, right.1, 1.0)
                .total_cmp(&claim(left.0, left.2, left.1, 1.0))
                .then(left.3.cmp(&right.3))
        });

        for (distance, losing, odds, contig) in shortlist {
            let Some(from) = of.get(&contig).copied().filter(|from| *from != id) else {
                continue;
            };
            let Some(donor) = bins.get(&from) else {
                continue;
            };
            if whole(quality.score(donor), settings) {
                continue;
            }
            if claim(distance, odds, losing, needed(quality, donor, contig)) < settings.confidence {
                continue;
            }
            let taken = with(&contigs, contig);
            if features.bin_size(&taken) > settings.max_bin_size {
                continue;
            }
            if quality.score(&taken).contamination > settings.contamination {
                continue;
            }

            let left = donor
                .iter()
                .copied()
                .filter(|other| *other != contig)
                .collect::<Vec<_>>();
            if left.is_empty() {
                bins.remove(&from);
                ledger.emptied += 1;
            } else {
                bins.insert(from, left);
            }
            of.insert(contig, id);
            ledger.taken += 1;
            ledger.taken_bp += features.length(contig);
            contigs = taken;
            moved += 1;
        }
        bins.insert(id, contigs);
    }
    moved
}

// The contigs that finish a cut genome carry no marker, so a marker delta cannot decide them
// and the move turns on which of the two bins the coverage overlap says owns the contig.
pub fn recruit(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    knn: &KnnGraph,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    settings: RecruitSettings,
) -> RecruitLedger {
    let metric = AggregateMetric::new(features.n_samples() * 2, features.distance_settings());
    let mut ledger = RecruitLedger::default();
    for _ in 0..settings.passes {
        ledger.passes += 1;
        let moved = pass(features, quality, knn, &metric, bins, settings, &mut ledger);
        if moved == 0 {
            break;
        }
    }
    ledger.bins = bins.len();
    ledger
}
