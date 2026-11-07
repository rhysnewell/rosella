#[path = "../support/bars.rs"]
mod bars;

use rosella::clustering::graph_partition::label_propagation;
use rosella::clustering::layout::base_level;
use rosella::clustering::leiden::leiden_from;
use rosella::embedding::Graph;
use rosella::embedding::features::ContigFeatures;
use rosella::markers::{ContigMarkers, MarkerSet};

fn leiden(graph: &Graph, gamma: f64, seed: u64) -> Vec<u32> {
    leiden_from(&base_level(graph, &vec![1.0; graph.rows()]), gamma, seed)
}

fn propagate(graph: &Graph, seed: u64) -> Vec<u32> {
    label_propagation(&base_level(graph, &vec![1.0; graph.rows()]), seed)
}

// Flat depth and no markers leave the arbiter's depth vote and eukaryote call nothing to act on.
pub(crate) struct Plain {
    coverage: ndarray::Array2<f64>,
    tnf: ndarray::Array2<f64>,
    lengths: Vec<usize>,
    pub(crate) markers: ContigMarkers,
}

impl Plain {
    pub(crate) fn new(contigs: usize) -> Self {
        Self {
            coverage: ndarray::Array2::zeros((contigs, 2)),
            tnf: ndarray::Array2::zeros((contigs, 2)),
            lengths: vec![1; contigs],
            markers: ContigMarkers::new(vec![Vec::new(); contigs], MarkerSet::parse_with("", None)),
        }
    }

    pub(crate) fn features(&self) -> ContigFeatures<'_> {
        ContigFeatures::new(&self.coverage, &self.tnf, &self.lengths)
    }
}

mod codelength_test;
mod combine_test;
mod conservation_test;
mod graph_partition_test;
mod leiden_test;
mod rung_bar_test;
