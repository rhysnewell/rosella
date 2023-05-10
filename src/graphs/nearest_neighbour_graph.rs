use annembed::fromhnsw::kgraph::KGraph;
use anyhow::Result;
use log::debug;
use num_traits::{FromPrimitive, Float};
use rayon::prelude::*;

/// Take a k-nearest neighbour graph and return a mutual k-nearest neighbour graph
/// A mutual k-nearest neighbour graph is a nearest neighbour graph where edges are only kept if they are mutual
/// i.e. if node A is a nearest neighbour of node B, and node B is a nearest neighbour of node A
pub fn mutual_knn<F: FromPrimitive + Float + std::fmt::UpperExp + Sync + Send + std::iter::Sum>(mut knn_graph: KGraph<F>) -> Result<KGraph<F>> {
    let mutual_nodes = knn_graph.get_neighbours()
        .into_par_iter()
        .enumerate()
        .map(|(node, neighbours)| {
            let mut mutual_neighbours = Vec::new();
            for neighbour in neighbours {
                if knn_graph.get_neighbours()[neighbour.node].iter().any(|edge| edge.node == node) {
                    mutual_neighbours.push(neighbour.clone());
                } else {
                    debug!("Node {} is a neighbour of node {}, but node {} is not a neighbour of node {}", neighbour.node, node, node, neighbour.node)
                }
            }
            Ok(mutual_neighbours)
        })
        .collect::<Result<Vec<_>>>()?;
    
    
    knn_graph.neighbours = mutual_nodes;
    Ok(knn_graph)
}


