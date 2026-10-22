#[path = "../support/bars.rs"]
mod bars;

use rosella::clustering::leiden::{base_level, leiden_from};
use rosella::embedding::Graph;

fn leiden(graph: &Graph, gamma: f64, seed: u64) -> Vec<i32> {
    leiden_from(&base_level(graph, None), gamma, seed)
}

mod codelength_test;
mod combine_test;
mod conservation_test;
mod graph_partition_test;
mod leiden_test;
mod rung_bar_test;
