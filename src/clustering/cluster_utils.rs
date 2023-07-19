use std::collections::HashMap;

use anyhow::Result;
use log::debug;
use ndarray::{ArrayBase, OwnedRepr, Dim, parallel::prelude::*, Axis};

pub fn condensed_pairwise_distance(data: &ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>) -> Vec<f64> {

    let distance_array = data
        .axis_iter(Axis(0))
        .into_par_iter()
        .enumerate().flat_map(|(i, x)| {
        data
            .axis_iter(Axis(0))
            .into_par_iter()
            .skip(i + 1)
            .map(|y| {
                // euclidean distance
                let mut d = x
                    .iter()
                    .zip(y.iter())
                    .map(|(x, y)| (x - y).powi(2))
                    .sum::<f64>()
                    .sqrt();
                if d.is_nan() {
                    d = 0.0;
                }
                d
            })
            .collect::<Vec<_>>()
    }).collect::<Vec<_>>();

    distance_array
}

/// Compute the silhoette score for this clustering result using the condensed distance matrix
/// returns the average silhouette score for all points, and optionally the silhouette score for each point
pub fn silhouette_score(distance_array: &[f64], cluster_labels: &HashMap<usize, Vec<usize>>, d: usize, score_each_point: bool) -> Result<(f64, Option<HashMap<usize, Vec<f64>>>)> {
    if cluster_labels.len() == 1 {
        return Ok((0.0, None));
    }

    let silhouette_scores: HashMap<usize, Vec<f64>> = cluster_labels
        .par_iter()
        .map(|(cluster_label, indices)| {
            let mut silhouettes = Vec::with_capacity(indices.len());
            // internal cluster distances
            for i in indices.iter() {
                // average distance between points in this cluster
                let mut a = 0.0;
    
                // minimum average distance from this cluster to antoher cluster
                let mut b = -1.0;
    
                let mut a_count = 0; // number of internal comparisons

                for j in indices.iter().skip(*i) {
                    if i == j {
                        continue;
                    }
                    let index = get_condensed_index(*i, *j, d);
                    a += distance_array[index];
                    a_count += 1;
                }

                if a_count == 0 {
                    a_count += 1;
                }
                a /= a_count as f64;
                // avearage distance to other clusters
                for (other_cluster_label, other_indices) in cluster_labels {
                    if other_cluster_label == cluster_label {
                        continue;
                    }
                    let mut b_sum = 0.0;
                    let mut b_count = 0; // number of external comparisons
                    
                    for j in other_indices {
                        let index = get_condensed_index(*i, *j, d);
                        b_sum += distance_array[index];
                        b_count += 1;
                    }
                    
    
                    if b_count == 0 {
                        b_count += 1;
                    }
    
                    b_sum /= b_count as f64;
                    if b_sum < b || b < 0.0 {
                        b = b_sum;
                    }
                }
                let s = (b - a) / a.max(b);
                if s.is_nan() {
                    debug!("a: {}, b: {}", a, b);
                    debug!("a_count: {}", a_count);
                }
    
                silhouettes.push(s);
            }
            (*cluster_label, silhouettes)
        })
        .collect::<HashMap<_, _>>();

    let score_sum = silhouette_scores
        .iter()
        .map(|(_, silhouettes)| {
            silhouettes.iter().sum::<f64>()
        })
        .sum::<f64>();

    let silhouette_scores = if score_each_point {
        Some(silhouette_scores)
    } else {
        None
    };

    Ok((score_sum / d as f64, silhouette_scores))
}


/// Get the condensed index for a given i, j, and d
/// i and j are the indices of the elements in the original matrix
/// d is the number of rows in the original matrix
pub fn get_condensed_index(i: usize, j: usize, d: usize) -> usize {
    if i == j {
        panic!("Cannot get condensed index for i == j");
    }

    let index = d * (d - 1) / 2 - (d - i) * (d - i - 1) / 2 + j - i - 1;
    index
}