use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet},
};

use anyhow::Result;
use hdbscan::{DistanceMetric, Hdbscan, HdbscanHyperParams, NnAlgorithm};
use log::{debug, trace};
use ndarray::{Array2, ArrayBase, Data, Ix2};
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::clustering::validity::dbcv;

/// flight sweeps min_cluster_size over ten values and keeps the best by validity. Its own
/// lower bound is computed but always collapses to 2, so the width is written out here.
const SWEEP_WIDTH: usize = 10;
const SMALLEST_CLUSTER: usize = 2;

/// Validity is quadratic in the points it scores, and the sweep scores every combination,
/// so it runs against a sample of a large embedding rather than all of it.
const VALIDITY_SAMPLE_LIMIT: usize = 5000;

pub struct HdbscanSettings {
    pub min_cluster_size: usize,
    pub min_samples: usize,
    pub seed: u64,
}

/// Cluster the embedding, sweeping the two size parameters and keeping the labelling with
/// the best density based cluster validity.
pub fn find_best_clusters<S: Data<Elem = f64> + Sync>(
    embeddings: &ArrayBase<S, Ix2>,
    seed: u64,
) -> Result<HDBSCANResult> {
    let rows = embeddings
        .rows()
        .into_iter()
        .map(|row| row.iter().map(|value| *value as f32).collect::<Vec<f32>>())
        .collect::<Vec<_>>();

    let sample = validity_sample(embeddings, seed);

    // The hdbscan crate reads the min_samples-th neighbour without checking there is one,
    // so a bin smaller than the sweep panics rather than erroring.
    let combinations = (SMALLEST_CLUSTER..SMALLEST_CLUSTER + SWEEP_WIDTH)
        .filter(|min_cluster_size| *min_cluster_size <= rows.len())
        .flat_map(|min_cluster_size| {
            (SMALLEST_CLUSTER..=min_cluster_size)
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
            let sampled_labels = sample
                .iter()
                .map(|index| labels[*index])
                .collect::<Vec<_>>();
            let validity = dbcv(&sample_rows(embeddings, &sample), &sampled_labels);

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

fn validity_sample<S: Data<Elem = f64>>(embeddings: &ArrayBase<S, Ix2>, seed: u64) -> Vec<usize> {
    let n = embeddings.nrows();
    if n <= VALIDITY_SAMPLE_LIMIT {
        return (0..n).collect();
    }

    let mut rng = StdRng::seed_from_u64(seed);
    let mut chosen = (0..n).collect::<Vec<_>>();
    for position in 0..VALIDITY_SAMPLE_LIMIT {
        chosen.swap(position, rng.random_range(position..n));
    }
    chosen.truncate(VALIDITY_SAMPLE_LIMIT);
    chosen.sort_unstable();
    chosen
}

fn sample_rows<S: Data<Elem = f64>>(
    embeddings: &ArrayBase<S, Ix2>,
    sample: &[usize],
) -> Array2<f64> {
    let mut rows = Array2::zeros((sample.len(), embeddings.ncols()));
    for (position, index) in sample.iter().enumerate() {
        rows.row_mut(position).assign(&embeddings.row(*index));
    }
    rows
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

    /// Fold another result in, renumbering its clusters so nothing collides.
    pub fn merge(&mut self, other: HDBSCANResult) {
        let mut next_cluster_id = self.cluster_map.keys().max().map_or(0, |id| id + 1);
        for indices in other.cluster_map.into_values() {
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
