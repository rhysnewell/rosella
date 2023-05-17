use std::{collections::HashMap, cmp::Ordering};

use anyhow::Result;
use env_logger::filter::Filter;
use log::{debug, info};
use ndarray::{ArrayBase, OwnedRepr, Dim};
use petal_clustering::{HDbscan, Fit};
use rayon::prelude::*;

use crate::clustering::cluster_utils::{condensed_pairwise_distance, silhouette_score};


/// Clusters the embeddings using HDBSCAN
pub fn find_best_clusters(embeddings: &ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>, starting_min_cluster_size: usize) -> Result<HDBSCANResult> {
    let n = embeddings.shape()[0]; // number of rows
    let end_size = starting_min_cluster_size + 10;

    let condensed_distances = condensed_pairwise_distance(embeddings);

    // for min_samples 5..15
    // calculate the best silhouette score, return the cluster map and outliers
    let mut cluster_results = (starting_min_cluster_size..end_size)
        .into_par_iter()
        .flat_map(|min_cluster_size| {
            (starting_min_cluster_size..end_size)
                .into_par_iter()
                .map(|min_samples| {
                    let mut clusterer = HDbscan::default();
                    clusterer.min_samples = min_samples;
                    clusterer.min_cluster_size = min_cluster_size;
                    clusterer.alpha = 1.0;
                    let (cluster_map, outliers) = clusterer.fit(&embeddings);
                    let (mut s_score, _) = silhouette_score(&condensed_distances, &cluster_map, n, false).expect("Failed to calculate silhouette score");
                    if s_score.is_nan() {
                        s_score = 0.0;
                    }
                    debug!("Min cluster size: {}, min samples: {}, silhouette score: {}", min_cluster_size, min_samples, s_score);
                    (cluster_map, outliers, s_score)
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    
    // sort by silhouette score, closer to 1 is better
    cluster_results.sort_by(|a, b| {
        a.2.partial_cmp(&b.2).unwrap_or(Ordering::Equal).reverse()
    });

    let (cluster_map, outliers, score) = cluster_results.remove(0);
    let (_, silhouette_scores) = silhouette_score(&condensed_distances, &cluster_map, n, true)?;
    info!("Best silhouette score: {}", score);
    let mut result = HDBSCANResult::new(cluster_map, outliers, score, silhouette_scores);

    result.find_deviant_points()?;
    // at this point the silhouette scores have been invalidated, due to the removal of points
    // so we need to recalculate them
    let (score, silhouette_scores) = silhouette_score(&condensed_distances, &result.cluster_map, n, true)?;
    result.score = score;
    result.silhouette_scores = silhouette_scores;

    Ok(result)
}

pub struct HDBSCANResult {
    pub cluster_map: HashMap<usize, Vec<usize>>,
    pub outliers: Vec<usize>,
    pub score: f64,
    pub silhouette_scores: Option<HashMap<usize, Vec<f64>>>,
}

impl HDBSCANResult {
    pub fn new(cluster_map: HashMap<usize, Vec<usize>>, outliers: Vec<usize>, score: f64, silhouette_scores: Option<HashMap<usize, Vec<f64>>>) -> Self {
        Self {
            cluster_map,
            outliers,
            score,
            silhouette_scores,
        }
    }

    /// Merge two HDBSCAN results together esnuring no points are duplicated
    /// and no cluster ids are duplicated. if there are duplicated cluster ids, update
    /// the new cluster ids to be unique
    /// additionally, merged the outlier contigs
    pub fn merge(&mut self, other: HDBSCANResult) {
        let mut minimum_cluster_id = self.cluster_map.keys().max().unwrap_or(&0) + 1;

        for (_, indices) in other.cluster_map.into_iter() {
            self.cluster_map.insert(minimum_cluster_id, indices);
            minimum_cluster_id += 1;
        }

        self.outliers = other.outliers;

        // silhouette scores are invalidated
        self.score = f64::NAN;
        self.silhouette_scores = None;
    }

    /// Reindex the clusters using a contig map. The map contains the current index
    /// as the key and the new index as value
    pub fn reindex_clusters(&mut self, contig_map: HashMap<usize, usize>) {
        let new_cluster_map = self.cluster_map.par_iter().map(|(cluster, points)| {
            let new_indices = points.par_iter().map(|point| {
                contig_map[point]
            }).collect::<Vec<usize>>();
            (*cluster, new_indices)
        }).collect::<HashMap<usize, Vec<usize>>>();

        self.cluster_map = new_cluster_map;

        let new_outliers = self.outliers.par_iter().map(|point| {
            contig_map[point]
        }).collect::<Vec<usize>>();

        self.outliers = new_outliers;
    }

    /// Find deviant points inside clusters by observing the silhouette score for all points
    /// and the silhouette score for inidividual points
    pub fn find_deviant_points(&mut self) -> Result<()> {
        let mut deviant_points = Vec::new();

        if let Some(silhouette_scores) = &self.silhouette_scores {
            // calculate the variance of silhouette scores
            // let n_clustered_points = self.cluster_map.values().map(|v| v.len()).sum::<usize>();
            // this represents the variance of all of the silhouette scores
            // kind of a like a base line background noise value
            // let variance = silhouette_scores.par_iter().map(|scores| {
            //     scores.iter().map(|score| {
            //         (score - self.score).powi(2)
            //     }).sum::<f64>()
            // }).sum::<f64>() / n_clustered_points as f64;


            let filtered_clusters = self.cluster_map.par_iter().filter_map(|(cluster, points)| {
                // calculate the variance of the silhouette scores for the cluster
                let mut cluster_variance = silhouette_scores[cluster].iter().map(|score| {
                    (score - self.score).powi(2)
                }).sum::<f64>() / points.len() as f64;
                cluster_variance = cluster_variance * 2.0;
                // this returns in the index in the point array, we need to switch this back to the actual values
                let (deviant_points_in_cluster, points_new): (Vec<_>, Vec<_>) = (0..points.len()).into_par_iter().partition(|point_index| {
                    let point_variance = (silhouette_scores[cluster][*point_index] - self.score).powi(2);
                    point_variance > cluster_variance
                });

                let points_new = points_new.into_par_iter().map(|point_index| {
                    points[point_index]
                }).collect::<Vec<usize>>();

                let deviant_points_in_cluster = deviant_points_in_cluster.into_par_iter().map(|point_index| {
                    points[point_index]
                }).collect::<Vec<usize>>();

                debug!("Cluster {} has {} deviant points variance {}", cluster, deviant_points_in_cluster.len(), cluster_variance);

                if deviant_points_in_cluster.len() > 0 {
                    Some(FilteredCluster {
                        cluster: *cluster,
                        points: points_new,
                        deviant_points: deviant_points_in_cluster,
                    })
                } else {
                    None
                }                
            }).collect::<Vec<FilteredCluster>>();

            // replace the cluster map with the filtered clusters and extend the deviant points
            filtered_clusters.into_iter().for_each(|filtered_cluster| {
                self.cluster_map.insert(filtered_cluster.cluster, filtered_cluster.points);
                deviant_points.par_extend(filtered_cluster.deviant_points);
            });

        } else {
            return Ok(())
        }


        self.outliers.par_extend(deviant_points);
        Ok(())
    }
}

struct FilteredCluster {
    cluster: usize,
    points: Vec<usize>,
    deviant_points: Vec<usize>,
}