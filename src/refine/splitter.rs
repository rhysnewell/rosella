use std::collections::{BTreeMap, HashMap, HashSet};

use log::debug;
use rayon::prelude::*;

use crate::{
    clustering::graph_partition::Partition,
    embedding::{features::ContigFeatures, knn::KnnSettings},
    quality::Scorer,
    refine::bar::describe_levels,
    refine::bin_stats::{BinStats, Thresholds, bin_stats},
    refine::gates::{Rejections, Trigger, TriggerCounts},
    refine::proposal::Proposal,
    refine::{bisect, cut_report::CutLog, floor},
    seeds::Seeds,
};

mod accept;

#[derive(Debug, Clone, Copy)]
pub struct RefineSettings {
    pub min_bin_size: usize,
    pub max_bin_size: usize,
    pub knn: KnnSettings,
    pub max_retries: usize,
    pub seeds: Seeds,
    pub max_contamination: Option<f64>,
    pub partition: Partition,
    pub split: bool,
    pub peel: bool,
}

pub struct Refined {
    pub bins: BTreeMap<usize, Vec<usize>>,
    pub unbinned: Vec<usize>,
    pub genome_floor: Option<usize>,
    pub cuts: Option<CutLog>,
    pub splits: usize,
}

#[derive(Default)]
struct Judges<'a> {
    assembly: Option<&'a crate::embedding::Graph>,
    quality: Option<&'a dyn Scorer>,
    contamination: HashMap<usize, f64>,
}

struct Round<'t> {
    thresholds: &'t Thresholds,
    eligible: usize,
}

#[derive(Default)]
struct Tally {
    rejections: Rejections,
    triggers: TriggerCounts,
}

pub struct Refiner<'a> {
    features: ContigFeatures<'a>,
    settings: RefineSettings,
    judges: Judges<'a>,
    refined: Refined,
    cached: BTreeMap<usize, BinStats>,
    survived: HashSet<usize>,
    next_bin_id: usize,
}

impl<'a> Refiner<'a> {
    pub fn new(
        features: ContigFeatures<'a>,
        settings: RefineSettings,
        bins: BTreeMap<usize, Vec<usize>>,
        unbinned: Vec<usize>,
    ) -> Self {
        let next_bin_id = bins.keys().max().map_or(1, |id| id + 1);
        Self {
            features,
            settings,
            judges: Judges::default(),
            refined: Refined {
                bins,
                unbinned,
                genome_floor: None,
                cuts: None,
                splits: 0,
            },
            cached: BTreeMap::new(),
            survived: HashSet::new(),
            next_bin_id,
        }
    }

    pub fn with_cuts(mut self, traced: bool) -> Self {
        self.refined.cuts = traced.then(CutLog::default);
        self
    }

    pub fn with_quality(mut self, quality: &'a dyn Scorer) -> Self {
        self.judges.quality = Some(quality);
        self
    }

    pub fn with_assembly(mut self, assembly: &'a crate::embedding::Graph) -> Self {
        self.judges.assembly = Some(assembly);
        self
    }

    pub fn with_contamination(mut self, contamination: HashMap<usize, f64>) -> Self {
        self.judges.contamination = contamination;
        self
    }

    fn measure_floor(&mut self) {
        self.refined.genome_floor = floor::floor(
            &self.features,
            &self.refined.bins,
            &self.refined.unbinned,
            self.settings.min_bin_size,
        );
        debug!(
            "Genome floor {:?}, split floor {}",
            self.refined.genome_floor,
            self.split_floor(),
        );
    }

    pub fn run(mut self) -> Refined {
        // Zero rounds means no refinement, so no result is reported for work that never ran.
        // The floor is still measured, or turning refinement off would move the pool's floor too.
        if self.settings.max_retries == 0 {
            self.measure_floor();
            return self.refined;
        }

        let mut total = Tally::default();
        for round in 0..self.settings.max_retries {
            let thresholds = self.refresh_stats();
            let pending = self
                .refined
                .bins
                .keys()
                .copied()
                .filter(|bin_id| !self.survived.contains(bin_id))
                .collect::<Vec<_>>();
            let eligible = pending
                .iter()
                .filter(|bin_id| {
                    bisect::eligible(
                        &self.features,
                        &self.refined.bins[bin_id],
                        self.settings.min_bin_size,
                    )
                })
                .count();
            self.measure_floor();
            let context = Round {
                thresholds: &thresholds,
                eligible,
            };
            // Proposing is the whole embed pipeline per bin and reads nothing another bin
            // writes, so it fans out. Applying stays in bin order, which is what fixes the ids.
            let progress = crate::progress::counted(
                crate::progress::Stage::RefiningBins,
                pending.len() as u64,
            );
            let proposals = pending
                .par_iter()
                .map(|bin_id| {
                    let proposal = self.propose(*bin_id, &context);
                    progress.inc(1);
                    proposal
                })
                .collect::<Vec<_>>();
            crate::progress::done(&progress);
            let mut tally = Tally::default();
            let mut split_this_round = 0;
            for (bin_id, proposal) in pending.into_iter().zip(proposals) {
                if self.apply(bin_id, proposal, &mut tally) {
                    split_this_round += 1;
                } else {
                    self.survived.insert(bin_id);
                }
            }

            debug!(
                "Refinement round {} split {} bins, turned away by {}",
                round, split_this_round, tally.rejections
            );
            debug!("Split levels are {}", describe_levels(&thresholds));
            debug!("Bins reached the bar as {}", tally.triggers);
            total.rejections.merge(&tally.rejections);
            total.triggers.merge(&tally.triggers);
            self.refined.splits += split_this_round;
            if split_this_round == 0 {
                break;
            }
        }
        debug!(
            "Refinement split {} bins, {} bins and {} unbinned contigs remain",
            self.refined.splits,
            self.refined.bins.len(),
            self.refined.unbinned.len()
        );
        debug!("Splits turned away by {}", total.rejections);
        debug!("Bins reached the bar as {}", total.triggers);
        self.refined
    }

    // Statistics for every bin, since the levels a bin is judged against are a quantile over
    // all of them. Bins that have not changed keep the figures they already had.
    fn refresh_stats(&mut self) -> Thresholds {
        let seed = self.settings.seeds.seed;
        let missing = self
            .refined
            .bins
            .iter()
            .filter(|(bin_id, _)| !self.cached.contains_key(bin_id))
            .map(|(bin_id, indices)| (*bin_id, indices.as_slice()))
            .collect::<Vec<_>>();

        {
            let _timer = crate::timing::scope("bin_stats");
            let fresh = missing
                .par_iter()
                .filter_map(|(bin_id, indices)| {
                    bin_stats(&self.features, indices, seed).map(|stats| (*bin_id, stats))
                })
                .collect::<Vec<_>>();
            self.cached.extend(fresh);
        }
        let bins = &self.refined.bins;
        self.cached.retain(|bin_id, _| bins.contains_key(bin_id));

        Thresholds::from_bins(self.cached.values())
    }

    fn apply(&mut self, bin_id: usize, proposal: Proposal, tally: &mut Tally) -> bool {
        let (trigger, outcome) = match proposal {
            Proposal::TooFewContigs => {
                tally.rejections.too_few_contigs += 1;
                return false;
            }
            Proposal::NoTrigger => {
                tally.rejections.no_trigger += 1;
                return false;
            }
            Proposal::NoClustering(trigger) => {
                tally.triggers.record(trigger);
                tally.rejections.no_clustering += 1;
                return false;
            }
            Proposal::Rejected(trigger, rejection) => {
                tally.triggers.record(trigger);
                tally.rejections.record(rejection);
                return false;
            }
            Proposal::Accepted(trigger, outcome) => {
                tally.triggers.record(trigger);
                (trigger, outcome)
            }
        };

        let members = self.refined.bins.remove(&bin_id).unwrap_or_default();
        self.cached.remove(&bin_id);
        if let Some(log) = self.refined.cuts.as_mut() {
            let features = &self.features;
            let pieces = crate::refine::owners::owners(outcome.kept.iter().enumerate());
            let stage = if trigger == Trigger::Peeled {
                "peel"
            } else {
                "split"
            };
            log.cut(stage, &members, &pieces, |contig| features.length(contig));
        }
        // A piece's figures were measured to accept the cut, so the next round reads them back.
        for (sub_bin, stats) in outcome.kept.into_iter().zip(outcome.stats) {
            if let Some(stats) = stats {
                self.cached.insert(self.next_bin_id, stats);
            }
            self.refined.bins.insert(self.next_bin_id, sub_bin);
            self.next_bin_id += 1;
        }
        self.refined.unbinned.extend(outcome.unbinned);
        true
    }
}
