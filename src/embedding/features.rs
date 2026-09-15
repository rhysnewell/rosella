use ndarray::Array2;

use crate::kmers::sketch::ContigSketches;
use crate::seeds::Seeds;

use crate::embedding::{
    Graph,
    knn::{KnnGraph, build_knn_with},
    metrics::{DistanceSettings, MIN_VAR, prepared::PreparedAggregate},
    fuzzy,
};

/// Coverage and composition for the whole assembly, addressed by contig index. Both the
/// initial embedding and refinement work through this, so neither owns the layout.
pub struct ContigFeatures<'a> {
    coverage: &'a Array2<f64>,
    tnf: &'a Array2<f64>,
    lengths: &'a [usize],
    distance: DistanceSettings,
    links: Option<&'a [(usize, usize)]>,
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

    /// An empty table is not the same as no table: it means the comparison ran and nothing
    /// aligned, which is the evidence that two genome-sized contigs are one genome.
    pub fn with_links(mut self, links: Option<&'a [(usize, usize)]>, weight: f32) -> Self {
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

    /// The variance floor `metabat` applies to this contig's coverage.
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

    /// Coverage and composition concatenated into one row per contig, which is what
    /// `AggregateMetric` splits back apart.
    pub fn rows(&self, indices: &[usize]) -> Vec<Vec<f64>> {
        indices
            .iter()
            .map(|index| {
                let coverage = self.coverage_row(*index);
                let tnf = self.tnf_row(*index);
                let mut row = Vec::with_capacity(coverage.len() + tnf.len());
                row.extend_from_slice(coverage);
                row.extend_from_slice(tnf);
                row
            })
            .collect()
    }

    pub(crate) fn floors(&self, indices: &[usize]) -> Vec<f64> {
        vec![MIN_VAR; indices.len()]
    }

    pub(crate) fn knn_size(&self, rows: usize, n_neighbours: usize) -> usize {
        let asked = if rows < n_neighbours * 10 {
            std::cmp::min(rows / 2, n_neighbours)
        } else {
            n_neighbours
        };
        asked.max(2)
    }

    fn combined_knn(
        &self,
        indices: &[usize],
        n_neighbours: usize,
        candidates: usize,
        seed: u64,
        stage: &'static str,
    ) -> KnnGraph {
        let _timer = crate::timing::scope(stage);
        let floors = self.floors(indices);
        let metric =
            PreparedAggregate::new(self.coverage, self.tnf, indices, &floors, self.distance);
        build_knn_with(
            indices.len(),
            self.knn_size(indices.len(), n_neighbours),
            candidates,
            seed,
            |a, b| metric.distance(a, b),
        )
    }

    pub fn contig_lengths(&self, indices: &[usize]) -> Vec<usize> {
        indices.iter().map(|index| self.lengths[*index]).collect()
    }

    pub fn graph_of(
        &self,
        indices: &[usize],
        n_neighbours: usize,
        seeds: Seeds,
        candidates: usize,
        stage: &'static str,
    ) -> Graph {
        let knn = self.knn_of(indices, n_neighbours, seeds, candidates, stage);
        self.graph_from_knn(indices, &knn)
    }

    pub fn knn_of(
        &self,
        indices: &[usize],
        n_neighbours: usize,
        seeds: Seeds,
        candidates: usize,
        stage: &'static str,
    ) -> KnnGraph {
        self.combined_knn(
            indices,
            n_neighbours,
            candidates,
            seeds.knn,
            stage,
        )
    }

    pub fn graph_from_knn(
        &self,
        indices: &[usize],
        knn: &KnnGraph,
    ) -> Graph {
        let graph = fuzzy::manifold_graph(indices.len(), knn, knn.indices.ncols());
        match self.links {
            Some(links) => crate::embedding::linked(graph, links, indices, self.link_weight),
            None => graph,
        }
    }
}


/// Rows of a standard-layout array are contiguous, so this never fails.
pub fn row_slice(array: &Array2<f64>, row: usize) -> &[f64] {
    array
        .row(row)
        .to_slice()
        .expect("array row is not contiguous")
}
