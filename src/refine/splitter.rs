use std::collections::{BTreeMap, HashMap, HashSet};

use log::{debug, info};
use rayon::prelude::*;

use crate::{
    clustering::{
        clusterer::{Partitioning, find_best_partition},
        objective::ClusterObjective,
    },
    embedding::features::ContigFeatures,
    refine::bar::{MIN_SPLIT_CONTIGS, describe_levels, should_split},
    refine::bin_stats::{AGGREGATE, BinStats, LevelSource, Thresholds, bin_stats},
    refine::gates::{Rejections, SplitRejection, Trigger, TriggerCounts},
    refine::proposal::{
        Proposal, SplitOutcome, contigs, judge_split, leaves_two_standing, tighter,
    },
    refine::{bisect, floor, peel},
};

const LEFTOVER_AGGREGATE: f64 = 0.5;

#[derive(Debug, Clone, Copy)]
pub struct RefineSettings {
    pub min_bin_size: usize,
    pub max_bin_size: usize,
    pub n_neighbours: usize,
    pub max_retries: usize,
    pub seeds: crate::seeds::Seeds,
    pub max_contamination: Option<f64>,
    pub overrides: crate::embedding::umap::EmbedOverrides,
    pub bisect: bool,
    pub levels: LevelSource,
    pub level_quantile: f64,
    pub partition: crate::clustering::graph_partition::Partition,
    pub node_size: crate::clustering::graph_partition::NodeSize,
    pub partition_resolution: Option<f64>,
    pub partition_theta: Option<f64>,
}

/// Splits chimeric bins by re-clustering them on their own.
pub struct Refiner<'a> {
    features: ContigFeatures<'a>,
    objective: &'a dyn ClusterObjective,
    settings: RefineSettings,
    assembly: Option<&'a crate::embedding::Graph>,
    pub bins: BTreeMap<usize, Vec<usize>>,
    pub unbinned: Vec<usize>,
    contamination: HashMap<usize, f64>,
    survived: HashSet<usize>,
    cached: BTreeMap<usize, BinStats>,
    next_bin_id: usize,
    eligible: usize,
    pub genome_floor: Option<usize>,
    rejections: Rejections,
    triggers: TriggerCounts,
}

impl<'a> Refiner<'a> {
    pub fn new(
        features: ContigFeatures<'a>,
        objective: &'a dyn ClusterObjective,
        settings: RefineSettings,
        bins: BTreeMap<usize, Vec<usize>>,
        unbinned: Vec<usize>,
    ) -> Self {
        let next_bin_id = bins.keys().max().map_or(1, |id| id + 1);
        Self {
            features,
            objective,
            settings,
            assembly: None,
            bins,
            unbinned,
            contamination: HashMap::new(),
            survived: HashSet::new(),
            cached: BTreeMap::new(),
            next_bin_id,
            eligible: 0,
            genome_floor: None,
            rejections: Rejections::default(),
            triggers: TriggerCounts::default(),
        }
    }

    pub fn with_assembly(mut self, assembly: &'a crate::embedding::Graph) -> Self {
        self.assembly = Some(assembly);
        self
    }

    pub fn with_contamination(mut self, contamination: HashMap<usize, f64>) -> Self {
        self.contamination = contamination;
        self
    }

    fn measure_floor(&mut self) {
        self.genome_floor = floor::floor(
            &self.features,
            &self.bins,
            &self.unbinned,
            self.settings.min_bin_size,
        );
    }

    pub fn run(&mut self) -> usize {
        // recover passes 0 rounds to mean no refinement, so falling through the loop would
        // report a refinement result for work that never ran. The floor is still measured,
        // or turning refinement off would silently move the pool's floor as well.
        if self.settings.max_retries == 0 {
            self.measure_floor();
            return 0;
        }

        let mut splits = 0;
        let mut rejections = Rejections::default();
        let mut triggers = TriggerCounts::default();
        for round in 0..self.settings.max_retries {
            // Both counters read per round, so a bin revisited across rounds was being counted
            // once per visit and the totals ran ahead of the bins they described.
            self.rejections = Rejections::default();
            self.triggers = TriggerCounts::default();
            let thresholds = self.refresh_stats();
            let mut split_this_round = 0;

            let pending = self
                .bins
                .keys()
                .copied()
                .filter(|bin_id| !self.survived.contains(bin_id))
                .collect::<Vec<_>>();
            self.eligible = pending
                .iter()
                .filter(|bin_id| {
                    bisect::eligible(&self.features, &self.bins[bin_id], self.settings.min_bin_size)
                })
                .count();
            self.measure_floor();
            // Proposing is the whole embed pipeline per bin and reads nothing another bin
            // writes, so it fans out. Applying stays in bin order, which is what fixes the ids.
            let proposals = pending
                .par_iter()
                .map(|bin_id| self.propose(*bin_id, &thresholds))
                .collect::<Vec<_>>();
            for (bin_id, proposal) in pending.into_iter().zip(proposals) {
                if self.apply(bin_id, proposal) {
                    split_this_round += 1;
                } else {
                    self.survived.insert(bin_id);
                }
            }

            debug!(
                "Refinement round {} split {} bins, turned away by {}",
                round, split_this_round, self.rejections
            );
            info!("Split levels are {}", describe_levels(&thresholds));
            debug!("Bins reached the bar as {}", self.triggers);
            rejections.merge(&self.rejections);
            triggers.merge(&self.triggers);
            splits += split_this_round;
            if split_this_round == 0 {
                break;
            }
        }
        self.rejections = rejections;
        self.triggers = triggers;

        info!(
            "Refinement split {} bins, {} bins and {} unbinned contigs remain",
            splits,
            self.bins.len(),
            self.unbinned.len()
        );
        info!("Splits turned away by {}", self.rejections);
        info!("Bins reached the bar as {}", self.triggers);
        splits
    }

    /// Statistics for every bin, since the levels a bin is judged against are an average
    /// over the large bins. Bins that have not changed keep the figures they already had.
    fn refresh_stats(&mut self) -> Thresholds {
        let seed = self.settings.seeds.sample;
        let missing = self
            .bins
            .iter()
            .filter(|(bin_id, _)| !self.cached.contains_key(bin_id))
            .map(|(bin_id, indices)| (*bin_id, indices.clone()))
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
        self.cached
            .retain(|bin_id, _| self.bins.contains_key(bin_id));

        Thresholds::from_bins(
            self.cached
                .iter()
                .map(|(bin_id, stats)| (self.features.bin_size(&self.bins[bin_id]), stats)),
            self.settings.levels,
            self.settings.level_quantile,
        )
    }

    fn propose(&self, bin_id: usize, thresholds: &Thresholds) -> Proposal {
        let indices = &self.bins[&bin_id];
        let Some(stats) = self.cached.get(&bin_id) else {
            return Proposal::NoStats;
        };
        let lengths = self.features.contig_lengths(indices);

        let partition = self.propose_partition(bin_id, indices, stats, &lengths, thresholds);
        if matches!(partition, Proposal::Accepted(..)) {
            return partition;
        }
        if let Some(outcome) = self.bisect(bin_id, indices) {
            return Proposal::Accepted(Trigger::Bisected, outcome);
        }
        match self.peel(indices, stats, &lengths) {
            Some(outcome) => {
                debug!(
                    "Peeled {} contigs out of bin {} of {}",
                    outcome.kept.len() - 1,
                    bin_id,
                    indices.len()
                );
                Proposal::Accepted(Trigger::Peeled, outcome)
            }
            None => partition,
        }
    }

    fn propose_partition(
        &self,
        bin_id: usize,
        indices: &[usize],
        stats: &BinStats,
        lengths: &[usize],
        thresholds: &Thresholds,
    ) -> Proposal {
        let bin_size = lengths.iter().sum::<usize>();
        let Some(trigger) = self.trigger(stats, lengths, bin_size, bin_id, thresholds)
        else {
            return if indices.len() < MIN_SPLIT_CONTIGS {
                Proposal::TooFewContigs
            } else {
                Proposal::NoTrigger
            };
        };
        self.attempt(bin_id, indices, stats, trigger)
    }

    fn attempt(
        &self,
        bin_id: usize,
        indices: &[usize],
        stats: &BinStats,
        trigger: Trigger,
    ) -> Proposal {
        let clustered = self.cluster_bin(indices, trigger == Trigger::Forced);
        let Some((result, validity)) = clustered else {
            return Proposal::NoClustering(trigger);
        };
        debug!(
            "Bin {} of {} contigs re-clustered into {} at validity {:.3}",
            bin_id,
            indices.len(),
            result.cluster_map.len(),
            validity
        );
        match self.accept(indices, stats, result) {
            Ok(outcome) => {
                debug!(
                    "Split bin {} of {} contigs into {}",
                    bin_id,
                    indices.len(),
                    outcome.kept.len()
                );
                Proposal::Accepted(trigger, outcome)
            }
            Err(rejection) => Proposal::Rejected(trigger, rejection),
        }
    }

    fn bisect(&self, bin_id: usize, indices: &[usize]) -> Option<SplitOutcome> {
        if !self.settings.bisect
            || !bisect::eligible(&self.features, indices, self.settings.min_bin_size)
        {
            return None;
        }
        let pieces = bisect::candidate(
            &self.features,
            indices,
            self.settings.min_bin_size,
            self.eligible,
            self.settings.seeds.sample,
        )?;
        debug!(
            "Bisected bin {} of {} contigs into {} and {}",
            bin_id,
            indices.len(),
            pieces[0].len(),
            pieces[1].len()
        );
        Some(SplitOutcome {
            kept: pieces.into(),
            unbinned: Vec::new(),
        })
    }

    /// The rest of the bin has to come out tighter once the lone contigs leave, weighted as
    /// if a contig on its own has no spread at all, which is what a genome in one contig is.
    fn peel(&self, indices: &[usize], stats: &BinStats, lengths: &[usize]) -> Option<SplitOutcome> {
        let peel = peel::candidate(indices, stats, lengths, self.settings.min_bin_size)?;
        let rest = bin_stats(&self.features, &peel.rest, self.settings.seeds.sample)?;
        let rest_bp = self.features.bin_size(&peel.rest) as f64;
        let lone_bp = self.features.bin_size(&peel.lone) as f64;
        if !tighter(
            rest.mean[AGGREGATE] * rest_bp / (rest_bp + lone_bp),
            stats,
            AGGREGATE,
        ) {
            return None;
        }
        let mut kept = peel
            .lone
            .into_iter()
            .map(|contig| vec![contig])
            .collect::<Vec<_>>();
        kept.push(peel.rest);
        Some(SplitOutcome {
            kept,
            unbinned: Vec::new(),
        })
    }

    fn apply(&mut self, bin_id: usize, proposal: Proposal) -> bool {
        let outcome = match proposal {
            Proposal::NoStats => {
                self.rejections.no_clustering += 1;
                return false;
            }
            Proposal::TooFewContigs => {
                self.rejections.too_few_contigs += 1;
                return false;
            }
            Proposal::NoTrigger => {
                self.rejections.no_trigger += 1;
                return false;
            }
            Proposal::NoClustering(trigger) => {
                self.triggers.record(trigger);
                self.rejections.no_clustering += 1;
                return false;
            }
            Proposal::Rejected(trigger, rejection) => {
                self.triggers.record(trigger);
                self.rejections.record(rejection);
                return false;
            }
            Proposal::Accepted(trigger, outcome) => {
                self.triggers.record(trigger);
                outcome
            }
        };

        self.bins.remove(&bin_id);
        self.cached.remove(&bin_id);
        for sub_bin in outcome.kept {
            self.bins.insert(self.next_bin_id, sub_bin);
            self.next_bin_id += 1;
        }
        self.unbinned.extend(outcome.unbinned);
        true
    }

    fn trigger(
        &self,
        stats: &BinStats,
        lengths: &[usize],
        bin_size: usize,
        bin_id: usize,
        thresholds: &Thresholds,
    ) -> Option<Trigger> {
        let over_budget = match (
            self.contamination.get(&bin_id),
            self.settings.max_contamination,
        ) {
            (Some(contamination), Some(budget)) => *contamination > budget,
            _ => false,
        };
        should_split(
            stats,
            lengths,
            bin_size,
            over_budget,
            self.settings.max_bin_size,
            thresholds,
        )
    }

    fn cluster_bin(&self, indices: &[usize], _forced: bool) -> Option<(Partitioning, f64)> {
        let seeds = self.settings.seeds;
        let graph = crate::refine::split_graph::bin_graph(
            &self.features,
            self.assembly,
            indices,
            self.settings.n_neighbours,
            seeds,
            &self.settings.overrides,
        );
        find_best_partition(
            &graph,
            &self.features.contig_lengths(indices),
            self.settings.node_size,
            self.objective,
            seeds.partition,
            self.settings.partition,
            self.settings.partition_resolution,
            self.settings.partition_theta,
        )
        .ok()
        .map(|result| {
            let validity = result.score;
            (result, validity)
        })
    }

    /// flight's `handle_new_embedding`, without its habit of leaving a rejected split's
    /// contigs in the unbinned list as well as in the bin they never left.
    fn accept(
        &self,
        indices: &[usize],
        stats: &BinStats,
        result: Partitioning,
    ) -> Result<SplitOutcome, SplitRejection> {
        let mut clusters = result
            .cluster_map
            .into_values()
            .map(|positions| contigs(indices, positions.into_iter()))
            .collect::<Vec<_>>();
        clusters.sort_unstable();
        let noise = contigs(indices, result.outliers.into_iter());

        let (kept, spare) = judge_split(clusters, noise, |cluster| {
            self.features.bin_size(cluster)
        })?;

        if !self.pieces_are_tighter(&kept, stats, AGGREGATE) {
            return Err(SplitRejection::NotTighter);
        }

        let outcome = self.place_leftovers(kept, spare);
        if !leaves_two_standing(&outcome.kept, self.split_floor(), |piece| {
                self.features.bin_size(piece)
            })
        {
            return Err(SplitRejection::Shredded);
        }
        if self.tests_modes()
            && !bisect::separates(
                &self.features,
                &outcome.kept,
                self.eligible,
                self.settings.seeds.sample,
            )
        {
            return Err(SplitRejection::Unimodal);
        }
        Ok(outcome)
    }

    /// A piece smaller than a genome is a shard of one, not a bin.
    fn split_floor(&self) -> usize {
        self.measured_genome().unwrap_or(self.settings.min_bin_size)
    }

    /// A floor under the bin floor is an estimate off one or two contigs, which says nothing
    /// about genome scale. `auto` reads that as a run with no closed genomes to measure.
    fn measured_genome(&self) -> Option<usize> {
        self.genome_floor
            .filter(|floor| *floor > self.settings.min_bin_size)
    }

    fn tests_modes(&self) -> bool {
        self.measured_genome().is_none()
    }

    /// Length weighted mean aggregate distance across the pieces against the whole. A
    /// chimeric bin falls apart into tighter pieces; a pure one does not.
    fn pieces_are_tighter(&self, kept: &[Vec<usize>], whole: &BinStats, column: usize) -> bool {
        let mut weighted = 0.0;
        let mut total = 0;
        for piece in kept.iter() {
            let Some(stats) = bin_stats(&self.features, piece, self.settings.seeds.sample) else {
                continue;
            };
            let size = self.features.bin_size(piece);
            weighted += stats.mean[column] * size as f64;
            total += size;
        }
        if total == 0 {
            return false;
        }

        tighter(weighted / total as f64, whole, column)
    }

    fn place_leftovers(&self, mut kept: Vec<Vec<usize>>, spare: Vec<usize>) -> SplitOutcome {
        let holds_together = bin_stats(&self.features, &spare, self.settings.seeds.sample)
            .is_some_and(|stats| stats.mean[AGGREGATE] <= LEFTOVER_AGGREGATE);

        if holds_together {
            kept.push(spare);
            return SplitOutcome {
                kept,
                unbinned: Vec::new(),
            };
        }

        SplitOutcome {
            kept,
            unbinned: spare,
        }
    }
}
