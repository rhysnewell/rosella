use std::collections::{BTreeMap, HashMap, HashSet};

use log::{debug, info, trace};
use ndarray::Array2;

use crate::{
    clustering::{
        clusterer::{HDBSCANResult, find_best_clusters},
        objective::{ClusterObjective, ScoreThresholds},
    },
    embedding::features::ContigFeatures,
    refine::bin_stats::{AGGREGATE, BinStats, EUCLIDEAN, METABAT, RHO, Thresholds, bin_stats},
};

/// Floors on each level, so a run where every bin looks alike does not start splitting on
/// noise. flight's `validate_bins`.
const FLOORS: [f64; 4] = [0.30, 0.15, 6.0, 0.35];
const MULTIPLIERS: [f64; 4] = [1.25, 1.5, 1.25, 1.5];

/// Contigs flagged as out of place have to add up to this before they alone trigger a
/// split.
const MISPLACED_LENGTH: usize = 1_000_000;

/// Leftovers this large become a bin of their own rather than going back to unbinned, so
/// long as they hold together.
const LEFTOVER_BIN_SIZE: usize = 200_000;
const LEFTOVER_AGGREGATE: f64 = 0.5;

/// Below this a bin cannot yield two sub-bins over `min_bin_size`, so the work is wasted.
const MIN_SPLIT_CONTIGS: usize = 10;

/// A bin whose mean coverage and composition distances average below this is already
/// clean, and is never looked at again.
const CLEAN_BIN: f64 = 0.05;

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
    pub seed: u64,
    pub max_contamination: Option<f64>,
    pub overrides: crate::embedding::umap::EmbedOverrides,
    pub largest_cluster: usize,
}

/// Splits chimeric bins by re-clustering them on their own. flight's `slow_refine`, minus
/// the KMeans fallback that only ever won because `validating.py:845` scored the HDBSCAN
/// branch off the wrong array.
/// The two bars a re-clustering has to clear. Both sit on the objective's scale, so they
/// move together when the objective changes.
#[derive(Debug, Clone, Copy)]
pub struct SplitBars {
    /// Derived per bin from how dirty it looks.
    pub target: f64,
    /// What a split into fewer than two clusters has to reach.
    pub single_cluster: f64,
}

pub struct Refiner<'a> {
    features: ContigFeatures<'a>,
    embedding: Option<&'a Array2<f64>>,
    objective: &'a dyn ClusterObjective,
    settings: RefineSettings,
    pub bins: BTreeMap<usize, Vec<usize>>,
    pub unbinned: Vec<usize>,
    contamination: HashMap<usize, f64>,
    survived: HashSet<usize>,
    cached: HashMap<usize, BinStats>,
    next_bin_id: usize,
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
            cached: HashMap::new(),
            next_bin_id,
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

            debug!("Refinement round {} split {} bins", round, split_this_round);
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
        splits
    }

    /// Statistics for every bin, since the levels a bin is judged against are an average
    /// over the large bins. Bins that have not changed keep the figures they already had.
    fn refresh_stats(&mut self) -> Thresholds {
        let seed = self.settings.seed;
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
            return false;
        };
        let stats = BinStats {
            mean: stats.mean,
            std: stats.std,
            per_contig: stats.per_contig.clone(),
        };
        let bin_size = self.features.bin_size(&indices);
        let Some(bars) = self.split_target(&stats, &indices, bin_size, bin_id, thresholds) else {
            return false;
        };

        let Some((result, validity)) = self.cluster_bin(&indices) else {
            return false;
        };
        trace!(
            "Bin {} of {} contigs re-clustered into {} at validity {:.3} against a bar of {:.3}",
            bin_id,
            indices.len(),
            result.cluster_map.len(),
            validity,
            bars.target
        );
        let Some(outcome) = self.accept(&indices, &stats, result, validity, bars) else {
            return false;
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
    ) -> Option<SplitBars> {
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

        let scale = self.objective.thresholds();
        min_validity(
            stats,
            &lengths,
            bin_size,
            over_budget,
            self.settings.max_bin_size,
            thresholds,
            scale,
        )
        .map(|target| SplitBars {
            target,
            single_cluster: scale.single_cluster,
        })
    }

    /// Cluster the bin where it already sits, then re-embed it on its own if that was not
    /// convincing. Whichever scores higher wins.
    fn cluster_bin(&self, indices: &[usize]) -> Option<(HDBSCANResult, f64)> {
        let seed = self.settings.seed;
        let mut best = self
            .embedding
            .map(|embedding| subset(embedding, indices))
            .and_then(|rows| {
                find_best_clusters(
                    &rows,
                    indices,
                    self.objective,
                    seed,
                    self.settings.largest_cluster,
                )
                .ok()
            })
            .map(|result| {
                let validity = result.score;
                (result, validity)
            });

        if best
            .as_ref()
            .is_none_or(|(_, validity)| *validity < self.objective.thresholds().re_embed_ceiling)
        {
            let embedded = self
                .features
                .embed(
                    indices,
                    self.settings.n_neighbours,
                    seed,
                    &self.settings.overrides,
                )
                .ok();
            let re_embedded = embedded
                .and_then(|embedded| {
                    find_best_clusters(
                        &embedded,
                        indices,
                        self.objective,
                        seed,
                        self.settings.largest_cluster,
                    )
                    .ok()
                })
                .map(|result| {
                    let validity = result.score;
                    (result, validity)
                });

            best = match (best, re_embedded) {
                (Some(first), Some(second)) if second.1 > first.1 => Some(second),
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
    ) -> Option<SplitOutcome> {
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
            self.settings.min_bin_size,
            |cluster| self.features.bin_size(cluster),
        )?;

        if !self.pieces_are_tighter(&kept, stats) {
            return None;
        }

        Some(self.place_leftovers(kept, spare))
    }

    /// Length weighted mean aggregate distance across the pieces against the whole. A
    /// chimeric bin falls apart into tighter pieces; a pure one does not.
    fn pieces_are_tighter(&self, kept: &[Vec<usize>], whole: &BinStats) -> bool {
        let mut weighted = 0.0;
        let mut total = 0;
        for piece in kept.iter() {
            let Some(stats) = bin_stats(&self.features, piece, self.settings.seed) else {
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
        let holds_together = self.features.bin_size(&spare) >= LEFTOVER_BIN_SIZE
            && bin_stats(&self.features, &spare, self.settings.seed)
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

/// The validity a re-clustering has to reach for a split to be taken, or `None` when the
/// bin looks fine as it stands. flight writes each level as
/// `min(max(floor, threshold * multiplier), mean + std * 1.5)`, but the second term always
/// exceeds the bin's own mean, so the trigger reduces to the first.
pub fn min_validity(
    stats: &BinStats,
    lengths: &[usize],
    bin_size: usize,
    over_budget: bool,
    max_bin_size: usize,
    thresholds: &Thresholds,
    scale: ScoreThresholds,
) -> Option<f64> {
    if lengths.len() < MIN_SPLIT_CONTIGS {
        return None;
    }
    if bin_size >= max_bin_size || over_budget {
        return Some(0.0);
    }

    let levels = levels(thresholds);
    let tripped = (0..4).any(|column| stats.mean[column] >= levels[column])
        || misplaced_length(stats, lengths, &levels) >= MISPLACED_LENGTH;

    let dirt = (stats.mean[METABAT] + stats.mean[RHO]) / 2.0;
    if tripped {
        return Some((1.0 - dirt).clamp(0.0, scale.tripped_ceiling));
    }
    if dirt > CLEAN_BIN {
        return Some((1.0 - dirt).clamp(0.0, scale.dirty_ceiling));
    }
    None
}

fn misplaced_length(stats: &BinStats, lengths: &[usize], levels: &[f64; 4]) -> usize {
    lengths
        .iter()
        .zip(stats.per_contig.iter())
        .filter(|(_, averages)| {
            averages[METABAT] >= levels[METABAT]
                || averages[RHO] >= levels[RHO]
                || averages[EUCLIDEAN] >= levels[EUCLIDEAN]
        })
        .map(|(length, _)| *length)
        .sum()
}

/// Which of a re-clustering's clusters are worth keeping, and what is left over. `None`
/// rejects the split outright and the original bin stands.
pub fn judge_split(
    clusters: Vec<Vec<usize>>,
    noise: Vec<usize>,
    validity: f64,
    bars: SplitBars,
    min_bin_size: usize,
    size_of: impl Fn(&[usize]) -> usize,
) -> Option<(Vec<Vec<usize>>, Vec<usize>)> {
    let distinct = clusters.len() + usize::from(!noise.is_empty());
    if distinct <= 1 || validity < bars.target {
        return None;
    }
    // Unreachable under a score that returns its worst value for a single cluster, which
    // DBCV does. A marker objective need not, so the guard stays.
    if clusters.len() < 2 && validity < bars.single_cluster {
        return None;
    }

    let bin_size = clusters
        .iter()
        .chain(std::iter::once(&noise))
        .map(|contigs| size_of(contigs))
        .sum::<usize>() as f64;
    if size_of(&noise) as f64 > MAX_NOISE_FRACTION * bin_size {
        return None;
    }

    let (kept, small): (Vec<_>, Vec<_>) = clusters
        .into_iter()
        .partition(|cluster| size_of(cluster) >= min_bin_size);
    if kept.is_empty() {
        return None;
    }

    let mut spare = noise;
    spare.extend(small.into_iter().flatten());
    spare.sort_unstable();
    Some((kept, spare))
}

fn levels(thresholds: &Thresholds) -> [f64; 4] {
    let mut levels = [0.0f64; 4];
    for column in 0..4 {
        levels[column] = FLOORS[column].max(MULTIPLIERS[column] * thresholds.mean[column]);
    }
    levels
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
