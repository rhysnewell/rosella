use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet},
};

use anyhow::Result;
use hdbscan::{DistanceMetric, Hdbscan, HdbscanHyperParams, NnAlgorithm};
use log::{debug, trace};
use ndarray::{ArrayBase, Data, Ix2};
use rayon::prelude::*;

use crate::clustering::objective::{ClusterObjective, EmbeddingSample};

/// flight sweeps min_cluster_size over ten values and keeps the best by validity. Its own
/// lower bound is computed but always collapses to 2, so the width is written out here.
pub const SWEEP_WIDTH: usize = 10;
const SMALLEST_CLUSTER: usize = 2;

/// The upper bound flight fixed for every assembly, 0.06% of a 60,000 contig assembly.
/// Raising it is measurably inert: DBCV's score falls monotonically from `min_cluster_size`
/// 3, so the larger values are tried and discarded. The bound is not what caps cluster size.
pub const DEFAULT_LARGEST_CLUSTER: usize = SMALLEST_CLUSTER + SWEEP_WIDTH - 1;

/// Cluster the embedding, sweeping the two size parameters and keeping the labelling the
/// objective scores highest. `contigs[i]` is the contig row `i` of `embeddings` came from.
pub fn find_best_clusters<S: Data<Elem = f64> + Sync>(
    embeddings: &ArrayBase<S, Ix2>,
    contigs: &[usize],
    objective: &dyn ClusterObjective,
    sample_seed: u64,
    largest_cluster: usize,
) -> Result<HDBSCANResult> {
    let _timer = crate::timing::scope("cluster");
    let rows = embeddings
        .rows()
        .into_iter()
        .map(|row| row.iter().map(|value| *value as f32).collect::<Vec<f32>>())
        .collect::<Vec<_>>();

    let sample = EmbeddingSample::new(embeddings.view(), sample_seed);

    // The hdbscan crate reads the min_samples-th neighbour without checking there is one,
    // so a bin smaller than the sweep panics rather than erroring.
    let combinations = cluster_sizes(largest_cluster)
        .into_iter()
        .filter(|min_cluster_size| *min_cluster_size <= rows.len())
        .flat_map(|min_cluster_size| {
            (SMALLEST_CLUSTER..=min_cluster_size.min(DEFAULT_LARGEST_CLUSTER))
                .map(move |min_samples| (min_cluster_size, min_samples))
        })
        .filter(|(_, min_samples)| *min_samples < rows.len())
        .collect::<Vec<_>>();

    let mut scored = combinations
        .par_iter()
        .filter_map(|(min_cluster_size, min_samples)| {
            let parameters = HdbscanHyperParams::builder()
                .min_cluster_size(*min_cluster_size)
                .min_samples(*min_samples)
                .dist_metric(DistanceMetric::Euclidean)
                .nn_algorithm(NnAlgorithm::Auto)
                .build();

            let labels = Hdbscan::new(&rows, parameters).cluster().ok()?;
            let validity = objective.score(&sample, contigs, &labels);

            trace!(
                "min_cluster_size {} min_samples {} validity {}",
                min_cluster_size, min_samples, validity
            );
            Some((labels, validity))
        })
        .collect::<Vec<_>>();

    if scored.is_empty() {
        anyhow::bail!("HDBSCAN failed for every parameter combination");
    }

    scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
    let (labels, validity) = scored.remove(0);
    debug!("Best validity {}", validity);

    Ok(HDBSCANResult::from_labels(&labels, validity))
}

/// The min_cluster_size values to try. Consecutive integers while they fit in the sweep
/// width, geometrically spaced beyond it, so raising the bound costs no extra HDBSCAN fits.
pub fn cluster_sizes(largest: usize) -> Vec<usize> {
    let largest = largest.max(SMALLEST_CLUSTER);
    if largest <= DEFAULT_LARGEST_CLUSTER {
        return (SMALLEST_CLUSTER..=largest).collect();
    }

    let ratio = largest as f64 / SMALLEST_CLUSTER as f64;
    let mut sizes = (0..SWEEP_WIDTH)
        .map(|step| {
            let fraction = step as f64 / (SWEEP_WIDTH - 1) as f64;
            (SMALLEST_CLUSTER as f64 * ratio.powf(fraction)).round() as usize
        })
        .collect::<Vec<_>>();
    sizes.dedup();
    sizes
}

pub struct HDBSCANResult {
    pub cluster_map: HashMap<usize, HashSet<usize>>,
    pub outliers: HashSet<usize>,
    pub score: f64,
}

impl HDBSCANResult {
    pub fn from_labels(labels: &[i32], score: f64) -> Self {
        let mut cluster_map: HashMap<usize, HashSet<usize>> = HashMap::new();
        let mut outliers = HashSet::new();

        for (index, label) in labels.iter().enumerate() {
            if *label < 0 {
                outliers.insert(index);
            } else {
                cluster_map
                    .entry(*label as usize)
                    .or_default()
                    .insert(index);
            }
        }

        Self {
            cluster_map,
            outliers,
            score,
        }
    }

    /// Fold another result in, renumbering its clusters so nothing collides. Ordered by
    /// lowest member, because hash order would give the same partition different bin names
    /// on every run.
    pub fn merge(&mut self, other: HDBSCANResult) {
        let mut next_cluster_id = self.cluster_map.keys().max().map_or(0, |id| id + 1);
        let mut incoming = other.cluster_map.into_values().collect::<Vec<_>>();
        incoming
            .sort_unstable_by_key(|indices| indices.iter().min().copied().unwrap_or(usize::MAX));
        for indices in incoming {
            self.cluster_map.insert(next_cluster_id, indices);
            next_cluster_id += 1;
        }
        self.outliers = other.outliers;
        self.score = f64::NAN;
    }

    /// Map positions within a subset back to their original contig indices.
    pub fn reindex_clusters(&mut self, contig_map: HashMap<usize, usize>) {
        self.cluster_map = self
            .cluster_map
            .par_iter()
            .map(|(cluster, points)| {
                let indices = points.par_iter().map(|point| contig_map[point]).collect();
                (*cluster, indices)
            })
            .collect();

        self.outliers = self
            .outliers
            .par_iter()
            .map(|point| contig_map[point])
            .collect();
    }
}
