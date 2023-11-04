use annembed::fromhnsw::kgraph::KGraph;
use anyhow::Result;
use log::debug;
use num_traits::{FromPrimitive, Float};
use rayon::prelude::*;

/// Take a k-nearest neighbour graph and return a mutual k-nearest neighbour graph
/// A mutual k-nearest neighbour graph is a nearest neighbour graph where edges are only kept if they are mutual
/// i.e. if node A is a nearest neighbour of node B, and node B is a nearest neighbour of node A
/// `keep_n_edges` is the number of edges to keep for each node, setting above 0 will ensure that each node has at least `keep_n_edges` edges
/// if they had `keep_n_edges` to begin with. Nodes that are kept are closest nodes as determined by the original k-nearest neighbour graph
pub fn mutual_knn<F: FromPrimitive + Float + std::fmt::UpperExp + Sync + Send + std::iter::Sum>(mut knn_graph: KGraph<F>, keep_n_edges: usize) -> Result<KGraph<F>> {
    let mutual_nodes = knn_graph.get_neighbours()
        .into_par_iter()
        .enumerate()
        .map(|(node, neighbours)| {
            let mut mutual_neighbours = Vec::new();
            for (neighbour_index, neighbour) in neighbours.iter().enumerate() {
                if neighbour_index < keep_n_edges {
                    mutual_neighbours.push(neighbour.clone());
                    continue;
                }

                if knn_graph.get_neighbours()[neighbour.node].iter().any(|edge| edge.node == node) {
                    mutual_neighbours.push(neighbour.clone());
                }
            }
            Ok(mutual_neighbours)
        })
        .collect::<Result<Vec<_>>>()?;
    
    
    // knn_graph.neighbours = mutual_nodes;
    Ok(knn_graph)
}

/// Intersect two graphs that were built using the same contig information
pub fn intersect<F: FromPrimitive + Float + std::fmt::UpperExp + Sync + Send + std::iter::Sum>(mut base_knn_graph: KGraph<F>, intersecting_graph: KGraph<F>, keep_n_edges: usize) -> Result<KGraph<F>> {

    // basically, iterate through the nodes in the base graph look at each nodes neighbours. Then ensure that the neighbour in the base graph is also a neighbour in the intersecting graph.
    // if not, then remove that neighbour from the base graph unless it is one of the keep_n_edges closest neighbours
    let intersecting_nodes = intersecting_graph.get_neighbours();
    let new_nodes = base_knn_graph.get_neighbours()
        .into_par_iter()
        .enumerate()
        .map(|(node, neighbours)| {
            let mut new_neighbours = Vec::new();
            let node_id_in_intersecting_graph = intersecting_graph.get_idx_from_dataid(base_knn_graph.get_data_id_from_idx(node).unwrap()).unwrap();
            for (neighbour_index, neighbour) in neighbours.iter().enumerate() {
                if neighbour_index < keep_n_edges {
                    new_neighbours.push(neighbour.clone());
                    continue;
                }

                match intersecting_graph.get_idx_from_dataid(base_knn_graph.get_data_id_from_idx(neighbour.node).unwrap()) {
                    Some(intersect_node_index) => {
                        if intersecting_nodes[intersect_node_index].iter().any(|edge| edge.node == node_id_in_intersecting_graph) {
                            new_neighbours.push(neighbour.clone());
                        }
                    },
                    None => {
                        debug!("Node {} in base graph has a neighbour {} that is not in the intersecting graph", node, neighbour.node);
                    }
                };
            }
            Ok(new_neighbours)
        })
        .collect::<Result<Vec<_>>>()?;

    // base_knn_graph.neighbours = new_nodes;
    Ok(base_knn_graph)
}


