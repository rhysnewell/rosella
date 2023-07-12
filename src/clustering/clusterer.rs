use std::{collections::HashMap, cmp::Ordering};

use annembed::fromhnsw::{kgraph::KGraph, kgraph_from_hnsw_all};
use anyhow::Result;
use hnsw_rs::prelude::{Hnsw, Distance};
use log::{debug, info, trace};
use ndarray::{ArrayBase, OwnedRepr, Dim, Dimension};
use petal_clustering::{HDbscan, Fit};
use rayon::prelude::*;

use crate::{clustering::cluster_utils::{condensed_pairwise_distance, silhouette_score}, embedding::embedder::{ContigInformation, DepthDistance}, graphs::nearest_neighbour_graph::mutual_knn};

#[derive(Debug, Clone)]
pub struct ClusterDistance;

impl<'a, D: Dimension> Distance<Vec<ContigInformation<'a, D>>> for ClusterDistance {
    fn eval(&self, va: &[Vec<ContigInformation<'a, D>>], vb: &[Vec<ContigInformation<'a, D>>]) -> f32 {
        let distance_sum = va.iter()
            .map(|contig_1| {
                vb.iter()
                    .map(|contig_2| {
                        let distance = DepthDistance.eval(contig_1, contig_2);
                        distance
                    })
                    .sum::<f32>()
            })
            .sum::<f32>();
        distance_sum / (va.len() * vb.len()) as f32
    }
}

/// Build a mutual K-NN graph, but keep at least `keep_n_edges` edges per node
/// Setting `keep_n_edges` to 0 can result in some nodes becoming disconnected, these nodes
/// should either be removed or reconnected to the graph
pub fn build_kgraph_of_clusters(cluster_information: &HashMap<usize, Vec<Vec<ContigInformation<'_, Dim<[usize; 1]>>>>>, keep_n_edges: usize, n_neighbours: usize, max_layers: usize, ef_construction: usize) -> Result<(KGraph<f64>, Vec<usize>)> {
    let mut contig_nn: Hnsw<Vec<ContigInformation<'_, Dim<[usize; 1]>>>, ClusterDistance> = Hnsw::new(
        n_neighbours, 
        cluster_information.len(), 
        max_layers, 
        ef_construction,
        ClusterDistance
    );
    contig_nn.set_keeping_pruned(true);

    info!("Inserting cluster data into HNSW.");
    let contig_data_for_insertion = (0..cluster_information.len())
        .into_iter()
        .map(|i| {
            let contig_information = cluster_information.get(&i).unwrap();
            (contig_information, i)
        })
        .collect::<Vec<_>>();
    debug!("Beginning insertion.");
    contig_nn.parallel_insert(&contig_data_for_insertion);
    
    info!("Constructing kgraph.");
    let mut kgraph: KGraph<f64> = kgraph_from_hnsw_all(&contig_nn, n_neighbours).unwrap();
    kgraph = mutual_knn(kgraph, keep_n_edges)?;
    let disconnected_nodes = kgraph
        .get_neighbours()
        .par_iter()
        .enumerate()
        .filter_map(|(node_index, nbrs)| {
            if nbrs.len() == 0 {
                debug!("Node {} has no neighbours.", node_index);
                Some(node_index)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    Ok((kgraph, disconnected_nodes))
}

/// Clusters the embeddings using HDBSCAN
pub fn find_best_clusters(embeddings: &ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>, starting_min_cluster_size: usize, starting_min_sample_size: usize) -> Result<HDBSCANResult> {
    let n = embeddings.shape()[0]; // number of rows
    let end_size = starting_min_cluster_size + 10;
    let end_sample_size = starting_min_sample_size + 10;

    let condensed_distances = condensed_pairwise_distance(embeddings);

    // for min_samples 5..15
    // calculate the best silhouette score, return the cluster map and outliers
    let mut cluster_results = (starting_min_cluster_size..end_size)
        .into_par_iter()
        .flat_map(|min_cluster_size| {
            (starting_min_sample_size..end_sample_size)
                .into_par_iter()
                .map(|min_samples| {
                    let mut clusterer = HDbscan::default();
                    clusterer.min_samples = min_samples;
                    clusterer.min_cluster_size = min_cluster_size;
                    clusterer.alpha = 1.0;
                    let (cluster_map, outliers) = clusterer.fit(&embeddings);
                    let cluster_map = renumber_clusters(cluster_map);
                    let (mut s_score, _) = silhouette_score(&condensed_distances, &cluster_map, n, false).expect("Failed to calculate silhouette score");
                    if s_score.is_nan() {
                        s_score = 0.0;
                    }
                    trace!("Min cluster size: {}, min samples: {}, silhouette score: {}", min_cluster_size, min_samples, s_score);
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
    trace!("Best silhouette score: {}", score);
    let result = HDBSCANResult::new(cluster_map, outliers, score, silhouette_scores);

    // result.find_deviant_points()?;
    // // at this point the silhouette scores have been invalidated, due to the removal of points
    // // so we need to recalculate them
    // let (score, silhouette_scores) = silhouette_score(&condensed_distances, &result.cluster_map, n, true)?;
    // result.score = score;
    // result.silhouette_scores = silhouette_scores;

    Ok(result)
}

/// Change the cluster ids from the original cluster ids to the new cluster ids that range from 0..n
/// where n is the number of clusters
pub fn renumber_clusters(cluster_map: HashMap<usize, Vec<usize>>) -> HashMap<usize, Vec<usize>> {
    let mut new_cluster_map = HashMap::new();
    let mut cluster_id = 0;
    for (_, points) in cluster_map.into_iter() {
        new_cluster_map.insert(cluster_id, points);
        cluster_id += 1;
    }
    new_cluster_map
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

    pub fn renumber_clusters(&mut self) {
        self.cluster_map = renumber_clusters(self.cluster_map.clone());
    }

    pub fn merge_cluster_from_map(&mut self, cluster_map: HashMap<usize, Vec<usize>>) {
        
        let mut current_cluster_id = self.cluster_map.keys().into_iter().max().unwrap() + 1;
        cluster_map.into_iter().for_each(|(_, indices)| {
            self.cluster_map.insert(current_cluster_id, indices);
            current_cluster_id += 1;
        });
    }

    /// We have clustered the clusters, now we take that result and merge any clusters that it says are conjoined
    pub fn merge_clusters_from_result(&mut self, clusters_of_clusters: Self) {
        let mut new_cluster_map = HashMap::with_capacity(self.cluster_map.len());

        debug!("N clusters {}", clusters_of_clusters.cluster_map.len());
        debug!("N outliers {}", clusters_of_clusters.outliers.len());

        let mut current_cluster_id = 0;
        clusters_of_clusters.cluster_map.into_iter().for_each(|(_, clusters_to_merge)| {
            let mut merged_indices = Vec::new();
            clusters_to_merge.into_iter().for_each(|cluster_to_merge| {
                let indices = self.cluster_map.remove(&cluster_to_merge).unwrap();
                merged_indices.extend(indices);
            });
            new_cluster_map.insert(current_cluster_id, merged_indices);
            current_cluster_id += 1;
        });

        // rescure the outliers too, as they are just clusters on their own
        clusters_of_clusters.outliers.into_iter().for_each(|outlier| {
            let indices = self.cluster_map.remove(&outlier).unwrap();
            new_cluster_map.insert(current_cluster_id, indices);
            current_cluster_id += 1;
        });


        self.cluster_map = new_cluster_map;
        self.score = clusters_of_clusters.score;
        self.silhouette_scores = None;
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

    /// This function prepares the HDBSCAN result for clustering of clusters
    /// it removes all clusters that are larger than the max_bin_size, and puts them aside to be reinserted later on
    /// All outliers are kept, unless they are less than the min_bin_size
    pub fn prepare_for_clustering_of_clusters(&mut self, max_bin_size: usize, min_bin_size: usize, contig_lengths: &[usize], remove_large_bins: bool) -> HashMap<usize, Vec<usize>> {
        let mut minimum_cluster_id = self.cluster_map.keys().max().unwrap_or(&0) + 1;

        // remove all clusters that are larger than the max_bin_size
        let clusters = std::mem::take(&mut self.cluster_map);
        let mut new_clusters = HashMap::with_capacity(clusters.len());
        for (cluster, points) in clusters.into_iter() {
            let bin_size = points.iter().map(|point| contig_lengths[*point]).sum::<usize>();
            if bin_size >= max_bin_size && remove_large_bins {
                // bin too big so we leave it out
                new_clusters.insert(cluster, points);
            } else {
                // bin is small enough so we reinsert it
                self.cluster_map.insert(cluster, points);
            }
        }

        let outliers = std::mem::take(&mut self.outliers);
        let mut new_outliers = Vec::with_capacity(outliers.len());
        for outlier in outliers.into_iter() {
            if contig_lengths[outlier] >= min_bin_size {
                self.cluster_map.insert(minimum_cluster_id, vec![outlier]);
                minimum_cluster_id += 1;
            } else {
                new_outliers.push(outlier);
            }
        }

        self.outliers = new_outliers;
        self.score = f64::NAN;
        self.silhouette_scores = None;

        return new_clusters;
    }

    /// Returns a HashMap with the contig id as key and the cluster id as value
    pub fn get_contig_to_cluster_map(&self) -> HashMap<usize, usize> {
        let clustered_contigs = self.cluster_map.par_iter().map(|(cluster, points)| {
            points.par_iter().map(|point| {
                (*point, *cluster)
            }).collect::<HashMap<usize, usize>>()
        }).flatten().collect::<HashMap<usize, usize>>();

        clustered_contigs
    }
}

struct FilteredCluster {
    cluster: usize,
    points: Vec<usize>,
    deviant_points: Vec<usize>,
}