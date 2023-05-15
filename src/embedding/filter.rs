use anyhow::Result;
use ndarray::{Array2, Data};
use linfa_nn::{CommonNearestNeighbour, distance::Distance, NearestNeighbour};
use hnsw_rs::hnsw::Hnsw;
use linfa::Float;

pub fn nearest_neighbour_filter(coverage_array: &Array2<f32>, kmer_array: &Array2<f32>) -> Vec<usize> {
    // CommonNearestNeighbour::KdTree
    Vec::new()
}

fn kd_tree_filter<F: Float, D: Distance<F> + Clone>(array: &Array2<F>, function: D, threshold: F) -> Result<Vec<usize>> {
    let mut disconnected_points = Vec::new();
    let nn = CommonNearestNeighbour::KdTree.from_batch(array, function.clone())?;
    for (i, row) in array.rows().into_iter().enumerate() {
        let nearest_neighbour = nn.k_nearest(row.view(), 1)?;
        // calculate the distance between the point and its nearest neighbour
        let distance = function.distance(nearest_neighbour[0].0.view(), row.view());
        if distance > threshold {
            disconnected_points.push(i);
        }
    }
    Ok(disconnected_points)
}