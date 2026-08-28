use anyhow::Result;
use ndarray::Array2;

use crate::embedding::{
    knn::{KnnGraph, build_knn},
    metrics::{AggregateMetric, DistanceSettings, aggregate_weight, variance_floor},
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
}

impl<'a> ContigFeatures<'a> {
    pub fn new(coverage: &'a Array2<f64>, tnf: &'a Array2<f64>, lengths: &'a [usize]) -> Self {
        Self {
            coverage,
            tnf,
            lengths,
            distance: DistanceSettings::default(),
            reference_length: median_length(lengths),
        }
    }

    pub fn with_distance(mut self, distance: DistanceSettings) -> Self {
        self.distance = distance;
        self
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

    pub fn weight(&self) -> f64 {
        aggregate_weight(self.n_samples())
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

    pub fn build_knn(
        &self,
        rows: &[Vec<f64>],
        indices: &[usize],
        n_neighbours: usize,
        seed: u64,
    ) -> KnnGraph {
        let _timer = crate::timing::scope("knn");
        let metric = AggregateMetric::new(self.coverage.ncols(), self.distance);
        let floors = indices
            .iter()
            .map(|index| self.variance_floor(*index))
            .collect::<Vec<_>>();
        let n_neighbours = if rows.len() < n_neighbours * 10 {
            std::cmp::min(rows.len() / 2, n_neighbours)
        } else {
            n_neighbours
        };
        build_knn(rows.len(), n_neighbours.max(2), seed, |a, b| {
            metric.distance(&rows[a], &rows[b], floors[a], floors[b])
        })
    }

    pub fn embed(
        &self,
        indices: &[usize],
        n_neighbours: usize,
        seed: u64,
        overrides: &umap::EmbedOverrides,
    ) -> Result<Array2<f64>> {
        let rows = self.rows(indices);
        let knn = self.build_knn(&rows, indices, n_neighbours, seed);
        let contig_lengths = indices
            .iter()
            .map(|index| self.lengths[*index])
            .collect::<Vec<_>>();

        let mut curve = umap::curve_params(&contig_lengths);
        curve.a = overrides.a.unwrap_or(curve.a);
        curve.b = overrides.b.unwrap_or(curve.b);

        let settings = umap::EmbedSettings {
            n_components: overrides
                .n_components
                .unwrap_or_else(|| umap::n_components(self.n_samples())),
            n_neighbours: knn.indices.ncols(),
            curve,
            n_epochs: umap::default_epochs(rows.len()),
            seed,
            vertex_weights: umap::length_weights(&contig_lengths, overrides.length_weight),
        };
        umap::embed(&rows, &knn, &settings)
    }
}

fn median_length(lengths: &[usize]) -> usize {
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
