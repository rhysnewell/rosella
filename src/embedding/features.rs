use anyhow::Result;
use log::info;
use ndarray::Array2;

use crate::seeds::Seeds;

use crate::embedding::{
    knn::{KnnGraph, build_knn},
    intersect,
    metrics::{AggregateMetric, DistanceSettings, View, ViewMetric, variance_floor},
    quality::neighbour_preservation,
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

    fn knn_with(
        &self,
        rows: &[Vec<f64>],
        indices: &[usize],
        n_neighbours: usize,
        seed: u64,
        distance: impl Fn(&[f64], &[f64], f64, f64) -> f64 + Sync,
    ) -> KnnGraph {
        let _timer = crate::timing::scope("knn");
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
            distance(&rows[a], &rows[b], floors[a], floors[b])
        })
    }

    fn combined_knn(
        &self,
        rows: &[Vec<f64>],
        indices: &[usize],
        n_neighbours: usize,
        seed: u64,
    ) -> KnnGraph {
        let metric = AggregateMetric::new(self.coverage.ncols(), self.distance);
        self.knn_with(rows, indices, n_neighbours, seed, move |a, b, x, y| {
            metric.distance(a, b, x, y)
        })
    }

    fn view_knn(
        &self,
        rows: &[Vec<f64>],
        indices: &[usize],
        view: View,
        n_neighbours: usize,
        seed: u64,
    ) -> KnnGraph {
        let metric = ViewMetric::new(
            self.coverage.ncols(),
            view,
            self.distance.aggregation,
            self.distance.presence_fraction,
        );
        self.knn_with(rows, indices, n_neighbours, seed, move |a, b, x, y| {
            metric.distance(a, b, x, y)
        })
    }

    pub fn embed(
        &self,
        indices: &[usize],
        n_neighbours: usize,
        seeds: Seeds,
        overrides: &umap::EmbedOverrides,
    ) -> Result<Array2<f64>> {
        let rows = self.rows(indices);
        let views = self.distance.views.selected();
        let graphs = if views.is_empty() {
            vec![self.combined_knn(&rows, indices, n_neighbours, seeds.knn)]
        } else {
            views
                .iter()
                .map(|view| self.view_knn(&rows, indices, *view, n_neighbours, seeds.knn))
                .collect()
        };

        let contig_lengths = indices
            .iter()
            .map(|index| self.lengths[*index])
            .collect::<Vec<_>>();

        // Every view has to be embedded into at least its own dimensionality, or the
        // spectral start is undetermined for the view that reads highest.
        let intrinsic_dimension = graphs
            .iter()
            .filter_map(|knn| knn.intrinsic_dimension())
            .fold(None, |widest: Option<f64>, estimate| {
                Some(widest.map_or(estimate, |value| value.max(estimate)))
            });

        let settings = umap::EmbedSettings {
            n_components: overrides
                .n_components
                .unwrap_or_else(|| umap::n_components(intrinsic_dimension, self.n_samples())),
            n_neighbours: graphs[0].indices.ncols(),
            curve: umap::Curve::from_overrides(&contig_lengths, overrides),
            n_epochs: overrides
                .n_epochs
                .unwrap_or_else(|| umap::default_epochs(rows.len())),
            seeds,
            vertex_weights: umap::length_weights(&contig_lengths, overrides.length_weight),
            spectral_init: overrides.spectral_init,
        };

        match intrinsic_dimension {
            Some(estimate) => info!(
                "Intrinsic dimensionality {estimate:.2} over {} contigs, embedding into {}",
                rows.len(),
                settings.n_components
            ),
            None => info!(
                "Intrinsic dimensionality not estimable over {} contigs, embedding into {}",
                rows.len(),
                settings.n_components
            ),
        }

        let mut manifolds = graphs
            .iter()
            .map(|knn| umap::manifold_graph(&rows, knn, &settings))
            .collect::<Vec<_>>();
        let curve = manifolds[0].1;
        let graph = if manifolds.len() == 1 {
            manifolds.remove(0).0
        } else {
            let _timer = crate::timing::scope("intersect");
            let learned = manifolds
                .into_iter()
                .map(|(graph, _)| graph)
                .collect::<Vec<_>>();
            intersect::intersect(&learned)
        };

        let embedding = umap::layout(&graph, curve, &settings)?;

        let reference = if views.is_empty() {
            graphs.into_iter().next()
        } else {
            Some(self.combined_knn(&rows, indices, n_neighbours, seeds.knn))
        };
        if let Some(kept) =
            reference.and_then(|knn| neighbour_preservation(&embedding, &knn, seeds.knn))
        {
            info!("Neighbour preservation {kept:.4} against the combined metric");
        }

        Ok(embedding)
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
