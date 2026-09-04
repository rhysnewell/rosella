use std::collections::{BTreeMap, HashMap, HashSet};

use log::{debug, info};
use ndarray::Array2;
use rayon::prelude::*;

use crate::{
    clustering::{
        clusterer::{HDBSCANResult, find_best_clusters, find_best_partition},
        objective::ClusterObjective,
    },
    embedding::features::ContigFeatures,
    refine::bar::{MIN_SPLIT_CONTIGS, describe_levels, min_validity},
    refine::bin_stats::{AGGREGATE, BinStats, LevelSource, Thresholds, bin_stats},
    refine::gates::{Rejections, SplitGate, SplitRejection, Trigger, TriggerCounts},
};

const LEFTOVER_AGGREGATE: f64 = 0.5;

/// At or under this the bin is already known to be bad, so the clustering it sits in is not
/// worth trusting however well it scores. flight's `validating.py:948`.
const REEMBED_BELOW_BAR: f64 = 0.5;

/// Under this the bin is being split whatever comes back, so the fresh embedding is taken
/// even when it scores lower. flight's `validating.py:1018`.
const FORCED_BAR: f64 = 0.1;

/// Noise above this fraction of the original bin means the split threw away more than it
/// explained.
const MAX_NOISE_FRACTION: f64 = 0.6;

/// The pieces have to be this much tighter than the bin they came out of. Density validity
/// says a labelling separates well, not that the bin was chimeric, and a pure genome
/// separates perfectly happily. Without this the split takes good bins apart.
const REQUIRED_IMPROVEMENT: f64 = 0.9;

#[derive(Debug, Clone, Copy)]
pub struct RefineSettings {
    pub min_bin_size: usize,
    pub max_bin_size: usize,
    pub n_neighbours: usize,
    pub max_retries: usize,
    pub seeds: crate::seeds::Seeds,
    pub max_contamination: Option<f64>,
    pub overrides: crate::embedding::umap::EmbedOverrides,
    pub largest_cluster: usize,
    pub gate: SplitGate,
    pub levels: LevelSource,
    pub level_quantile: f64,
    pub split_bar: Option<f64>,
    pub partition: crate::clustering::graph_partition::Partition,
    pub partition_resolution: Option<f64>,
    pub partition_theta: Option<f64>,
}

/// `single_cluster` is on the objective's scale. `target` is in distance units and is
/// compared against a validity anyway, which is flight's conflation, not the port's.
#[derive(Debug, Clone, Copy)]
pub struct SplitBars {
    /// Derived per bin from how dirty it looks.
    pub target: f64,
    /// What a split into fewer than two clusters has to reach.
    pub single_cluster: f64,
}

/// Splits chimeric bins by re-clustering them on their own. flight's `slow_refine`, minus
/// the KMeans fallback that only ever won because `validating.py:845` scored the HDBSCAN
/// branch off the wrong array.
pub struct Refiner<'a> {
    features: ContigFeatures<'a>,
    embedding: Option<&'a Array2<f64>>,
    objective: &'a dyn ClusterObjective,
    settings: RefineSettings,
    pub bins: BTreeMap<usize, Vec<usize>>,
    pub unbinned: Vec<usize>,
    contamination: HashMap<usize, f64>,
    survived: HashSet<usize>,
    cached: BTreeMap<usize, BinStats>,
    next_bin_id: usize,
    rejections: Rejections,
    triggers: TriggerCounts,
}

impl<'a> Refiner<'a> {
    pub fn new(
        features: ContigFeatures<'a>,
        embedding: Option<&'a Array2<f64>>,
        objective: &'a dyn ClusterObjective,
        settings: RefineSettings,
        bins: BTreeMap<usize, Vec<usize>>,
        unbinned: Vec<usize>,
    ) -> Self {
        let next_bin_id = bins.keys().max().map_or(1, |id| id + 1);
        Self {
            features,
            embedding,
            objective,
            settings,
            bins,
            unbinned,
            contamination: HashMap::new(),
            survived: HashSet::new(),
            cached: BTreeMap::new(),
            next_bin_id,
            rejections: Rejections::default(),
            triggers: TriggerCounts::default(),
        }
    }

    pub fn with_contamination(mut self, contamination: HashMap<usize, f64>) -> Self {
        self.contamination = contamination;
        self
    }

    pub fn run(&mut self) -> usize {
        // recover passes 0 rounds to mean no refinement, so falling through the loop would
        // report a refinement result for work that never ran.
        if self.settings.max_retries == 0 {
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
        let bin_size = self.features.bin_size(indices);
        let Some((bars, trigger)) = self.split_target(stats, indices, bin_size, bin_id, thresholds)
        else {
            return if indices.len() < MIN_SPLIT_CONTIGS {
                Proposal::TooFewContigs
            } else {
                Proposal::NoTrigger
            };
        };

        let Some((result, validity)) = self.cluster_bin(indices, bars.target) else {
            return Proposal::NoClustering(trigger);
        };
        debug!(
            "Bin {} of {} contigs re-clustered into {} at validity {:.3} against a bar of {:.3}",
            bin_id,
            indices.len(),
            result.cluster_map.len(),
            validity,
            bars.target
        );
        match self.accept(indices, stats, result, validity, bars) {
            Ok(outcome) => {
                debug!(
                    "Split bin {} of {} contigs into {} at validity {:.3} against {:.3}",
                    bin_id,
                    indices.len(),
                    outcome.kept.len(),
                    validity,
                    bars.target
                );
                Proposal::Accepted(trigger, outcome)
            }
            Err(rejection) => Proposal::Rejected(trigger, rejection),
        }
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

    fn split_target(
        &self,
        stats: &BinStats,
        indices: &[usize],
        bin_size: usize,
        bin_id: usize,
        thresholds: &Thresholds,
    ) -> Option<(SplitBars, Trigger)> {
        let over_budget = match (
            self.contamination.get(&bin_id),
            self.settings.max_contamination,
        ) {
            (Some(contamination), Some(budget)) => *contamination > budget,
            _ => false,
        };
        let lengths = indices
            .iter()
            .map(|index| self.features.length(*index))
            .collect::<Vec<_>>();

        min_validity(
            stats,
            &lengths,
            bin_size,
            over_budget,
            self.settings.max_bin_size,
            thresholds,
        )
        .map(|(target, trigger)| {
            let bars = SplitBars {
                target: match (self.settings.split_bar, trigger) {
                    (Some(bar), Trigger::Tripped { .. }) => bar,
                    _ => target,
                },
                single_cluster: self.objective.thresholds().single_cluster,
            };
            (bars, trigger)
        })
    }

    /// Cluster the bin where it already sits, then re-embed it on its own if that was not
    /// convincing. Whichever scores higher wins, and a tie goes to the fresh embedding.
    fn cluster_bin(&self, indices: &[usize], target: f64) -> Option<(HDBSCANResult, f64)> {
        let seeds = self.settings.seeds;
        let scored = |result: HDBSCANResult| {
            let validity = result.score;
            (result, validity)
        };
        let mut best = self
            .embedding
            .map(|embedding| subset(embedding, indices))
            .and_then(|rows| {
                find_best_clusters(
                    &rows,
                    indices,
                    self.objective,
                    seeds.sample,
                    self.settings.largest_cluster,
                )
                .ok()
            })
            .map(scored);

        let unconvincing = best
            .as_ref()
            .is_none_or(|(_, validity)| *validity < self.objective.thresholds().re_embed_ceiling);

        if unconvincing || target <= REEMBED_BELOW_BAR {
            let graph_only =
                self.settings.partition.reads_graph() && !self.objective.needs_layout();
            let embedded = if graph_only {
                let graph = self.features.graph_of(
                    indices,
                    self.settings.n_neighbours,
                    seeds,
                    &self.settings.overrides,
                );
                Some((None, graph))
            } else {
                self.features
                    .embed_with_graph(
                        indices,
                        self.settings.n_neighbours,
                        seeds,
                        &self.settings.overrides,
                    )
                    .ok()
                    .map(|(embedded, graph)| (Some(embedded), graph))
            };
            let re_embedded = embedded
                .and_then(|(embedded, graph)| {
                    if self.settings.partition.reads_graph() {
                        find_best_partition(
                            &graph,
                            embedded.as_ref(),
                            indices,
                            self.objective,
                            seeds.sample,
                            seeds.partition,
                            self.settings.partition,
                            self.settings.partition_resolution,
                            self.settings.partition_theta,
                        )
                        .ok()
                    } else {
                        find_best_clusters(
                            embedded.as_ref()?,
                            indices,
                            self.objective,
                            seeds.sample,
                            self.settings.largest_cluster,
                        )
                        .ok()
                    }
                })
                .map(scored);

            best = match (best, re_embedded) {
                (Some(first), Some(second)) if second.1 >= first.1 => Some(second),
                (Some(first), Some(second))
                    if target < FORCED_BAR && first.1 <= REEMBED_BELOW_BAR =>
                {
                    Some(second)
                }
                (Some(first), _) => Some(first),
                (None, second) => second,
            };
        }

        best
    }

    /// flight's `handle_new_embedding`, without its habit of leaving a rejected split's
    /// contigs in the unbinned list as well as in the bin they never left.
    fn accept(
        &self,
        indices: &[usize],
        stats: &BinStats,
        result: HDBSCANResult,
        validity: f64,
        bars: SplitBars,
    ) -> Result<SplitOutcome, SplitRejection> {
        let mut clusters = result
            .cluster_map
            .into_values()
            .map(|positions| contigs(indices, positions.into_iter()))
            .collect::<Vec<_>>();
        clusters.sort_unstable();
        let noise = contigs(indices, result.outliers.into_iter());

        let (kept, spare) = judge_split(
            clusters,
            noise,
            validity,
            bars,
            self.settings.gate,
            |cluster| self.features.bin_size(cluster),
        )?;

        if self.settings.gate.is_strict() && !self.pieces_are_tighter(&kept, stats) {
            return Err(SplitRejection::NotTighter);
        }

        Ok(self.place_leftovers(kept, spare))
    }

    /// Length weighted mean aggregate distance across the pieces against the whole. A
    /// chimeric bin falls apart into tighter pieces; a pure one does not.
    fn pieces_are_tighter(&self, kept: &[Vec<usize>], whole: &BinStats) -> bool {
        let mut weighted = 0.0;
        let mut total = 0;
        for piece in kept.iter() {
            let Some(stats) = bin_stats(&self.features, piece, self.settings.seeds.sample) else {
                continue;
            };
            let size = self.features.bin_size(piece);
            weighted += stats.mean[AGGREGATE] * size as f64;
            total += size;
        }
        if total == 0 {
            return false;
        }

        weighted / total as f64 <= whole.mean[AGGREGATE] * REQUIRED_IMPROVEMENT
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

enum Proposal {
    NoStats,
    TooFewContigs,
    NoTrigger,
    NoClustering(Trigger),
    Rejected(Trigger, SplitRejection),
    Accepted(Trigger, SplitOutcome),
}

struct SplitOutcome {
    kept: Vec<Vec<usize>>,
    unbinned: Vec<usize>,
}

/// Size is not a bar here. A piece too small to write out can still recruit or merge its way
/// over the floor, so `bin_writer` applies `min_bin_size` once, at the end.
pub fn judge_split(
    clusters: Vec<Vec<usize>>,
    noise: Vec<usize>,
    validity: f64,
    bars: SplitBars,
    gate: SplitGate,
    size_of: impl Fn(&[usize]) -> usize,
) -> Result<(Vec<Vec<usize>>, Vec<usize>), SplitRejection> {
    let distinct = clusters.len() + usize::from(!noise.is_empty());
    if distinct <= 1 {
        return Err(SplitRejection::SingleCluster);
    }
    if validity < bars.target {
        return Err(SplitRejection::BelowTarget);
    }
    // Unreachable under a score that returns its worst value for a single cluster, which
    // DBCV does. A marker objective need not, so the guard stays.
    if clusters.len() < 2 && validity < bars.single_cluster {
        return Err(SplitRejection::SingleCluster);
    }

    let bin_size = clusters
        .iter()
        .chain(std::iter::once(&noise))
        .map(|contigs| size_of(contigs))
        .sum::<usize>() as f64;
    if gate.is_strict() && size_of(&noise) as f64 > MAX_NOISE_FRACTION * bin_size {
        return Err(SplitRejection::AllNoise);
    }

    Ok((clusters, noise))
}

fn contigs(indices: &[usize], positions: impl Iterator<Item = usize>) -> Vec<usize> {
    let mut contigs = positions
        .map(|position| indices[position])
        .collect::<Vec<_>>();
    contigs.sort_unstable();
    contigs
}

fn subset(embedding: &Array2<f64>, indices: &[usize]) -> Array2<f64> {
    let mut rows = Array2::zeros((indices.len(), embedding.ncols()));
    for (position, index) in indices.iter().enumerate() {
        rows.row_mut(position).assign(&embedding.row(*index));
    }
    rows
}
