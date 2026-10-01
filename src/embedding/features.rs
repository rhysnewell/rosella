use ndarray::Array2;

use crate::kmers::sketch::ContigSketches;

use crate::embedding::{
    Graph, fuzzy,
    knn::{KnnGraph, build_knn_with},
    metrics::{DistanceSettings, Point, prepared::PreparedAggregate},
};

// The initial embedding and refinement both work through this, so neither owns the layout.
pub struct ContigFeatures<'a> {
    coverage: &'a Array2<f64>,
    tnf: &'a Array2<f64>,
    lengths: &'a [usize],
    distance: DistanceSettings,
    links: Option<&'a [crate::assembly_graph::Link]>,
    link_weight: f32,
    sketches: Option<&'a ContigSketches>,
}

impl<'a> ContigFeatures<'a> {
    pub fn new(coverage: &'a Array2<f64>, tnf: &'a Array2<f64>, lengths: &'a [usize]) -> Self {
        Self {
            coverage,
            tnf,
            lengths,
            distance: DistanceSettings::default(),
            links: None,
            link_weight: 0.0,
            sketches: None,
        }
    }

    pub fn with_distance(mut self, distance: DistanceSettings) -> Self {
        self.distance = distance;
        self
    }

    pub fn with_links(
        mut self,
        links: Option<&'a [crate::assembly_graph::Link]>,
        weight: f32,
    ) -> Self {
        self.links = links.filter(|_| weight > 0.0);
        self.link_weight = weight;
        self
    }

    pub fn with_sketches(mut self, sketches: Option<&'a ContigSketches>) -> Self {
        self.sketches = sketches;
        self
    }

    pub fn sketches(&self) -> Option<&'a ContigSketches> {
        self.sketches
    }

    pub fn distance_settings(&self) -> DistanceSettings {
        self.distance
    }

    pub fn n_samples(&self) -> usize {
        self.coverage.ncols() / 2
    }

    pub fn coverage_row(&self, index: usize) -> &[f64] {
        row_slice(self.coverage, index)
    }

    pub fn tnf_row(&self, index: usize) -> &[f64] {
        row_slice(self.tnf, index)
    }

    pub fn length(&self, index: usize) -> usize {
        self.lengths[index]
    }

    pub fn bin_size(&self, indices: &[usize]) -> usize {
        indices.iter().map(|index| self.lengths[*index]).sum()
    }

    pub fn point(&self, index: usize) -> Point {
        Point::new(
            self.coverage_row(index),
            self.tnf_row(index),
            self.distance.presence_fraction,
        )
    }

    pub fn points(&self, indices: &[usize]) -> Vec<Point> {
        indices.iter().map(|index| self.point(*index)).collect()
    }

    pub(crate) fn knn_size(&self, rows: usize, n_neighbours: usize) -> usize {
        let asked = if rows < n_neighbours * 10 {
            std::cmp::min(rows / 2, n_neighbours)
        } else {
            n_neighbours
        };
        asked.max(2)
    }

    pub fn prepared(&self, indices: &[usize]) -> PreparedAggregate {
        PreparedAggregate::new(self.coverage, self.tnf, indices, self.distance)
    }

    pub fn knn_of(
        &self,
        indices: &[usize],
        n_neighbours: usize,
        candidates: usize,
        seed: u64,
        stage: &'static str,
    ) -> KnnGraph {
        let _timer = crate::timing::scope(stage);
        let metric = self.prepared(indices);
        build_knn_with(
            indices.len(),
            self.knn_size(indices.len(), n_neighbours),
            candidates,
            seed,
            &metric,
        )
    }

    pub fn contig_lengths(&self, indices: &[usize]) -> Vec<usize> {
        indices.iter().map(|index| self.lengths[*index]).collect()
    }

    pub fn graph_of(
        &self,
        indices: &[usize],
        n_neighbours: usize,
        candidates: usize,
        seed: u64,
        stage: &'static str,
    ) -> Graph {
        let knn = self.knn_of(indices, n_neighbours, candidates, seed, stage);
        self.graph_from_knn(indices, &knn)
    }

    pub fn graph_from_knn(&self, indices: &[usize], knn: &KnnGraph) -> Graph {
        let graph = fuzzy::manifold_graph(indices.len(), knn);
        self.linked(graph, indices)
    }

    fn linked(&self, graph: Graph, indices: &[usize]) -> Graph {
        match self.links {
            Some(links) => crate::embedding::linked(graph, links, indices, self.link_weight),
            None => graph,
        }
    }
}

// Rows of a standard-layout array are contiguous, so this never fails.
pub fn row_slice(array: &Array2<f64>, row: usize) -> &[f64] {
    array
        .row(row)
        .to_slice()
        .expect("array row is not contiguous")
}
