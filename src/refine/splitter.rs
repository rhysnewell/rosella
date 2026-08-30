use std::collections::{BTreeMap, HashMap, HashSet};

use log::{debug, info};
use ndarray::Array2;

use crate::{
    clustering::{
        clusterer::{HDBSCANResult, find_best_clusters},
        objective::ClusterObjective,
    },
    embedding::features::ContigFeatures,
    refine::bar::{MIN_SPLIT_CONTIGS, min_validity},
    refine::bin_stats::{AGGREGATE, BinStats, Thresholds, bin_stats},
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
        for round in 0..self.settings.max_retries {
            let thresholds = self.refresh_stats();
            let mut split_this_round = 0;

            for bin_id in self.bins.keys().copied().collect::<Vec<_>>() {
                if self.survived.contains(&bin_id) {
                    continue;
                }
                if self.visit(bin_id, &thresholds) {
                    split_this_round += 1;
                } else {
                    self.survived.insert(bin_id);
                }
            }

            debug!(
                "Refinement round {} split {} bins, turned away by {}",
                round, split_this_round, self.rejections
            );
            debug!("Bins reached the bar as {}", self.triggers);
            splits += split_this_round;
            if split_this_round == 0 {
                break;
            }
        }

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
            for (bin_id, indices) in missing {
                if let Some(stats) = bin_stats(&self.features, &indices, seed) {
                    self.cached.insert(bin_id, stats);
                }
            }
        }
        self.cached
            .retain(|bin_id, _| self.bins.contains_key(bin_id));

        Thresholds::from_bins(
            self.cached
                .iter()
                .map(|(bin_id, stats)| (self.features.bin_size(&self.bins[bin_id]), stats)),
        )
    }

    fn visit(&mut self, bin_id: usize, thresholds: &Thresholds) -> bool {
        let indices = self.bins[&bin_id].clone();
        let Some(stats) = self.cached.get(&bin_id) else {
            self.rejections.no_clustering += 1;
            return false;
        };
        let stats = BinStats {
            mean: stats.mean,
            std: stats.std,
            per_contig: stats.per_contig.clone(),
        };
        let bin_size = self.features.bin_size(&indices);
        let Some((bars, trigger)) = self.split_target(&stats, &indices, bin_size, bin_id, thresholds)
        else {
            if indices.len() < MIN_SPLIT_CONTIGS {
                self.rejections.too_few_contigs += 1;
            } else {
                self.rejections.no_trigger += 1;
            }
            return false;
        };
        self.triggers.record(trigger);

        let Some((result, validity)) = self.cluster_bin(&indices, bars.target) else {
            self.rejections.no_clustering += 1;
            return false;
        };
        debug!(
            "Bin {} of {} contigs re-clustered into {} at validity {:.3} against a bar of {:.3}",
            bin_id,
            indices.len(),
            result.cluster_map.len(),
            validity,
            bars.target
        );
        let outcome = match self.accept(&indices, &stats, result, validity, bars) {
            Ok(outcome) => outcome,
            Err(rejection) => {
                self.rejections.record(rejection);
                return false;
            }
        };

        debug!(
            "Split bin {} of {} contigs into {} at validity {:.3} against {:.3}",
            bin_id,
            indices.len(),
            outcome.kept.len(),
            validity,
            bars.target
        );

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
                target,
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
            let embedded = self
                .features
                .embed(
                    indices,
                    self.settings.n_neighbours,
                    seeds,
                    &self.settings.overrides,
                )
                .ok();
            let re_embedded = embedded
                .and_then(|embedded| {
                    find_best_clusters(
                        &embedded,
                        indices,
                        self.objective,
                        seeds.sample,
                        self.settings.largest_cluster,
                    )
                    .ok()
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
