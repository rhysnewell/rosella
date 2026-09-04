use anyhow::Result;
use log::info;
use ndarray::Array2;

use crate::seeds::Seeds;

use crate::embedding::{
    Graph, intersect,
    knn::{KnnGraph, MAX_CANDIDATES, build_knn, build_knn_with},
    manifold::{self, GraphWeights},
    metrics::{DistanceSettings, View, ViewMetric, prepared::PreparedAggregate, variance_floor},
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
        let metric = PreparedAggregate::new(rows, &floors, self.coverage.ncols(), self.distance);
        build_knn_with(
            rows.len(),
            self.knn_size(rows.len(), n_neighbours),
            candidates.unwrap_or(MAX_CANDIDATES),
            seed,
            |a, b| metric.distance(a, b),
        )
    }

    fn view_knn(
        &self,
        rows: &[Vec<f64>],
        indices: &[usize],
        view: View,
        n_neighbours: usize,
        seed: u64,
    ) -> KnnGraph {
        let _timer = crate::timing::scope("knn");
        let floors = self.floors(indices);
        let metric = ViewMetric::new(
            self.coverage.ncols(),
            view,
            self.distance.aggregation,
            self.distance.presence_fraction,
        );
        build_knn(
            rows.len(),
            self.knn_size(rows.len(), n_neighbours),
            seed,
            |a, b| metric.distance(&rows[a], &rows[b], floors[a], floors[b]),
        )
    }

    pub fn embed(
        &self,
        indices: &[usize],
        n_neighbours: usize,
        seeds: Seeds,
        overrides: &umap::EmbedOverrides,
    ) -> Result<Array2<f64>> {
        self.embed_with_graph(indices, n_neighbours, seeds, overrides)
            .map(|(embedding, _)| embedding)
    }

    pub fn embed_with_graph(
        &self,
        indices: &[usize],
        n_neighbours: usize,
        seeds: Seeds,
        overrides: &umap::EmbedOverrides,
    ) -> Result<(Array2<f64>, Graph)> {
        let manifold = self.manifold_of(indices, n_neighbours, seeds, overrides);
        let contig_lengths = self.contig_lengths(indices);

        // Every view has to be embedded into at least its own dimensionality, or the
        // spectral start is undetermined for the view that reads highest.
        let intrinsic_dimension = manifold
            .knn
            .iter()
            .filter_map(|knn| knn.intrinsic_dimension())
            .fold(None, |widest: Option<f64>, estimate| {
                Some(widest.map_or(estimate, |value| value.max(estimate)))
            });

        let settings = umap::EmbedSettings {
            n_components: overrides
                .n_components
                .unwrap_or_else(|| umap::n_components(intrinsic_dimension, self.n_samples())),
            n_epochs: overrides
                .n_epochs
                .unwrap_or_else(|| umap::default_epochs(indices.len())),
            seeds,
            vertex_weights: umap::length_weights(&contig_lengths, overrides.length_weight),
            spectral_init: overrides.spectral_init,
        };

        match intrinsic_dimension {
            Some(estimate) => info!(
                "Intrinsic dimensionality {estimate:.2} over {} contigs, embedding into {}",
                indices.len(),
                settings.n_components
            ),
            None => info!(
                "Intrinsic dimensionality not estimable over {} contigs, embedding into {}",
                indices.len(),
                settings.n_components
            ),
        }

        let embedding = umap::layout(&manifold.graph, manifold.curve, &settings)?;

        if let Some(kept) = self
            .preservation_knn(indices, n_neighbours, overrides, seeds, manifold.knn)
            .and_then(|knn| neighbour_preservation(&embedding, &knn, seeds.knn))
        {
            info!("Neighbour preservation {kept:.4} against the combined metric");
        }

        Ok((embedding, manifold.graph))
    }

    fn preservation_knn(
        &self,
        indices: &[usize],
        n_neighbours: usize,
        overrides: &umap::EmbedOverrides,
        seeds: Seeds,
        built: Vec<KnnGraph>,
    ) -> Option<KnnGraph> {
        if !overrides.report_preservation {
            return None;
        }
        if self.distance.views.selected().is_empty() {
            return built.into_iter().next();
        }
        Some(self.combined_knn(
            &self.rows(indices),
            indices,
            n_neighbours,
            overrides.knn_candidates,
            seeds.knn,
        ))
    }

    fn contig_lengths(&self, indices: &[usize]) -> Vec<usize> {
        indices.iter().map(|index| self.lengths[*index]).collect()
    }

    /// A graph partition ranked on a graph score reads the manifold and nothing else, and the
    /// layout on top of it is the most expensive stage in the run.
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
        let views = self.distance.views.selected();
        if views.is_empty() {
            vec![self.combined_knn(
                rows,
                indices,
                n_neighbours,
                overrides.knn_candidates,
                seeds.knn,
            )]
        } else {
            views
                .iter()
                .map(|view| self.view_knn(rows, indices, *view, n_neighbours, seeds.knn))
                .collect()
        }
    }

    fn manifold_of(
        &self,
        indices: &[usize],
        n_neighbours: usize,
        seeds: Seeds,
        overrides: &umap::EmbedOverrides,
    ) -> Manifold {
        let rows = self.rows(indices);
        let knn = self.knn_views(&rows, indices, n_neighbours, seeds, overrides);

        let contig_lengths = self.contig_lengths(indices);
        let curve = umap::Curve::from_overrides(&contig_lengths, overrides);
        let pinned = match curve {
            umap::Curve::Pinned(curve) => curve,
            umap::Curve::Fit { .. } => umap::curve_params(&contig_lengths),
        };
        let width = knn[0].indices.ncols();
        let mut manifolds = knn
            .iter()
            .map(|knn| match overrides.graph_weights {
                GraphWeights::Fuzzy => umap::manifold_graph(rows.len(), knn, width, curve),
                GraphWeights::Snn => (manifold::shared_neighbours(knn), pinned),
                GraphWeights::LocalScale => (manifold::local_scaled(knn), pinned),
            })
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

        Manifold { graph, curve, knn }
    }
}

struct Manifold {
    graph: Graph,
    curve: umap::CurveParams,
    knn: Vec<KnnGraph>,
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
