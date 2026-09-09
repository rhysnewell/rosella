use ndarray::Array2;

use crate::kmers::sketch::ContigSketches;
use crate::seeds::Seeds;

use crate::embedding::{
    Graph,
    knn::{KnnGraph, MAX_CANDIDATES, build_knn_with},
    manifold::{self, GraphWeights},
    metrics::{
        CompositionMetric, DistanceSettings, euclidean,
        prepared::PreparedAggregate, variance_floor,
    },
    umap,
};

/// Coverage and composition for the whole assembly, addressed by contig index. Both the
/// initial embedding and refinement work through this, so neither owns the layout.
pub struct ContigFeatures<'a> {
    coverage: &'a Array2<f64>,
    tnf: &'a Array2<f64>,
    lengths: &'a [usize],
    distance: DistanceSettings,
    reference_length: usize,
    sketches: Option<&'a ContigSketches>,
}

impl<'a> ContigFeatures<'a> {
    pub fn new(coverage: &'a Array2<f64>, tnf: &'a Array2<f64>, lengths: &'a [usize]) -> Self {
        Self {
            coverage,
            tnf,
            lengths,
            distance: DistanceSettings::default(),
            reference_length: median_length(lengths),
            sketches: None,
        }
    }

    pub fn with_distance(mut self, distance: DistanceSettings) -> Self {
        self.distance = distance;
        self.distance.composition_scale = composition_scale(self.tnf, distance.composition);
        self
    }

    /// An empty table is not the same as no table: it means the comparison ran and nothing
    /// aligned, which is the evidence that two genome-sized contigs are one genome.
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
    pub fn variance_floor(&self, index: usize) -> f64 {
        variance_floor(
            self.lengths[index],
            self.reference_length,
            self.distance.length_scaled_variance,
        )
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

    fn floors(&self, indices: &[usize]) -> Vec<f64> {
        indices
            .iter()
            .map(|index| self.variance_floor(*index))
            .collect()
    }

    fn knn_size(&self, rows: usize, n_neighbours: usize) -> usize {
        let asked = if rows < n_neighbours * 10 {
            std::cmp::min(rows / 2, n_neighbours)
        } else {
            n_neighbours
        };
        asked.max(2)
    }

    fn combined_knn(
        &self,
        rows: &[Vec<f64>],
        indices: &[usize],
        n_neighbours: usize,
        candidates: Option<usize>,
        seed: u64,
    ) -> KnnGraph {
        let _timer = crate::timing::scope("knn");
        let floors = self.floors(indices);
        let metric = PreparedAggregate::new(
            rows,
            &floors,
            self.coverage.ncols(),
            self.distance,
        );
        build_knn_with(
            rows.len(),
            self.knn_size(rows.len(), n_neighbours),
            candidates.unwrap_or(MAX_CANDIDATES),
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
        overrides: &umap::EmbedOverrides,
    ) -> Graph {
        self.manifold_of(indices, n_neighbours, seeds, overrides)
            .graph
    }

    pub fn knn_of(
        &self,
        indices: &[usize],
        n_neighbours: usize,
        seeds: Seeds,
        overrides: &umap::EmbedOverrides,
    ) -> Vec<KnnGraph> {
        self.knn_views(&self.rows(indices), indices, n_neighbours, seeds, overrides)
    }

    fn knn_views(
        &self,
        rows: &[Vec<f64>],
        indices: &[usize],
        n_neighbours: usize,
        seeds: Seeds,
        overrides: &umap::EmbedOverrides,
    ) -> Vec<KnnGraph> {
        vec![self.combined_knn(
            rows,
            indices,
            n_neighbours,
            overrides.knn_candidates,
            seeds.knn,
        )]
    }

    /// Read mapping puts strain siblings at one depth, so a bin the coverage view calls uniform
    /// has only composition left to split it on.
    fn manifold_of(
        &self,
        indices: &[usize],
        n_neighbours: usize,
        seeds: Seeds,
        overrides: &umap::EmbedOverrides,
    ) -> Manifold {
        let rows = self.rows(indices);
        let knn = self.knn_views(&rows, indices, n_neighbours, seeds, overrides);
        self.manifold_from(indices, &rows, knn, overrides)
    }

    fn manifold_from(
        &self,
        indices: &[usize],
        rows: &[Vec<f64>],
        knn: Vec<KnnGraph>,
        overrides: &umap::EmbedOverrides,
    ) -> Manifold {
        let contig_lengths = self.contig_lengths(indices);
        let curve = umap::Curve::from_overrides(&contig_lengths, overrides);
        let _pinned = match curve {
            umap::Curve::Pinned(curve) => curve,
            umap::Curve::Fit { .. } => umap::curve_params(&contig_lengths),
        };
        let width = knn[0].indices.ncols();
        let graph = match overrides.graph_weights {
            GraphWeights::Fuzzy => umap::manifold_graph(rows.len(), &knn[0], width, curve).0,
            GraphWeights::Snn => manifold::shared_neighbours(&knn[0]),
            GraphWeights::LocalScale => manifold::local_scaled(&knn[0]),
        };

        Manifold { graph }
    }
}

struct Manifold {
    graph: Graph,
}

/// Aitchison distance has no natural ceiling, and the refiner's bars are calibrated to rho's
/// [0, 2]. Mapping the run's own median pair to 1 puts the two on the same footing.
const SCALE_ROWS: usize = 256;

fn composition_scale(tnf: &Array2<f64>, metric: CompositionMetric) -> f64 {
    if metric != CompositionMetric::Aitchison || tnf.nrows() < 2 {
        return 1.0;
    }
    let sample = (0..tnf.nrows())
        .step_by(tnf.nrows().div_ceil(SCALE_ROWS))
        .collect::<Vec<_>>();
    let mut distances = Vec::with_capacity(sample.len() * sample.len() / 2);
    for (position, a) in sample.iter().enumerate() {
        for b in &sample[position + 1..] {
            distances.push(euclidean(row_slice(tnf, *a), row_slice(tnf, *b)));
        }
    }
    distances.sort_by(f64::total_cmp);
    match distances.get(distances.len() / 2) {
        Some(median) if *median > 0.0 => *median,
        _ => 1.0,
    }
}

pub fn median_length(lengths: &[usize]) -> usize {
    if lengths.is_empty() {
        return 1;
    }
    let mut sorted = lengths.to_vec();
    sorted.sort_unstable();
    sorted[sorted.len() / 2].max(1)
}

/// Rows of a standard-layout array are contiguous, so this never fails.
pub fn row_slice(array: &Array2<f64>, row: usize) -> &[f64] {
    array
        .row(row)
        .to_slice()
        .expect("array row is not contiguous")
}
