#[path = "../support/bars.rs"]
mod bars;

use rosella::clustering::graph_partition::label_propagation;
use rosella::clustering::layout::base_level;
use rosella::clustering::leiden::leiden_from;
use rosella::embedding::Graph;

fn leiden(graph: &Graph, gamma: f64, seed: u64) -> Vec<u32> {
    leiden_from(&base_level(graph, &vec![1.0; graph.rows()]), gamma, seed)
}

fn propagate(graph: &Graph, seed: u64) -> Vec<u32> {
    label_propagation(&base_level(graph, &vec![1.0; graph.rows()]), seed)
}

mod codelength_test;
mod combine_test;
mod conservation_test;
mod graph_partition_test;
mod leiden_test;
mod rung_bar_test;
