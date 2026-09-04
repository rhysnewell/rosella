mod search;
mod state;

use crate::clustering::graph_partition::{compact, label_propagation, node_degrees};
use crate::clustering::leiden::Level;
use crate::embedding::Graph;

/// Degree corrected block model under minimum description length (Peixoto, 2014). Seeded from
/// label propagation, because a move costs the blocks a node reaches and singletons maximise that.
pub fn sbm(graph: &Graph, seed: u64) -> Vec<i32> {
    let nodes = graph.rows();
    if nodes == 0 {
        return Vec::new();
    }

    let degrees = node_degrees(graph);
    if degrees.iter().sum::<f64>() <= 0.0 {
        return vec![0; nodes];
    }

    let start = label_propagation(graph, seed);
    let block_count = start
        .iter()
        .copied()
        .max()
        .map_or(1, |top| top as usize + 1);
    let of = start
        .iter()
        .map(|label| *label as usize)
        .collect::<Vec<_>>();

    let level = Level::from_graph(graph);
    let blocks = search::search(&level.neighbours, &degrees, of, block_count, seed);

    compact(&blocks.iter().map(|block| *block as i32).collect::<Vec<_>>())
}
