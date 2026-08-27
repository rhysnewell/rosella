use anyhow::Result;
use log::debug;
use ndarray::Array2;
use umap_rs::{
    EuclideanMetric, GraphParams, ManifoldParams, MetricType, OptimizationParams, Optimizer, Umap,
    UmapConfig,
};

use crate::embedding::{knn::KnnGraph, spectral::spectral_init};

const MIN_COMPONENTS: usize = 2;
const SMALL_DATASET: usize = 10_000;
const MAX_COMPONENTS: usize = 10;

/// Parameters of UMAP's distance to probability curve, `1 / (1 + a * x^(2b))`.
#[derive(Debug, Clone, Copy)]
pub struct CurveParams {
    pub a: f32,
    pub b: f32,
}

/// Contig length at which the cumulative length of the shortest contigs first passes
/// `percent` of the assembly. Matches flight's `nX`, which sorts ascending.
pub fn n_x(contig_lengths: &[usize], percent: f64) -> usize {
    if contig_lengths.is_empty() {
        return 0;
    }

    let mut lengths = contig_lengths.to_vec();
    lengths.sort_unstable();

    let target = lengths.iter().sum::<usize>() as f64 * (percent / 100.0);
    let mut running = 0.0;
    for length in lengths.iter() {
        running += *length as f64;
        if running > target {
            return *length;
        }
    }
    *lengths.last().unwrap()
}

/// flight derives the curve from assembly contiguity rather than exposing it, so that a
/// fragmented assembly gets a looser embedding than a contiguous one.
pub fn curve_params(contig_lengths: &[usize]) -> CurveParams {
    let log10 = |value: usize| (value.max(1) as f64).log10();

    let numerator = log10(n_x(contig_lengths, 25.0)).max(50_000f64.log10());
    let denominator = 500_000f64.log10().max(log10(n_x(contig_lengths, 75.0)));
    let b = (0.1 * (numerator / denominator) + 0.2).clamp(0.3, 0.4);
    let a = (log10(n_x(contig_lengths, 10.0)) * 0.1 + 1.0).clamp(1.4, 2.0);

    CurveParams {
        a: a as f32,
        b: b as f32,
    }
}

pub fn n_components(n_samples: usize) -> usize {
    n_samples.clamp(MIN_COMPONENTS, MAX_COMPONENTS)
}

pub struct EmbedSettings {
    pub n_components: usize,
    pub n_neighbours: usize,
    pub curve: CurveParams,
    pub n_epochs: usize,
    pub seed: u64,
}

pub fn default_epochs(n_points: usize) -> usize {
    if n_points <= SMALL_DATASET { 500 } else { 200 }
}

pub fn embed(rows: &[Vec<f64>], knn: &KnnGraph, settings: &EmbedSettings) -> Result<Array2<f64>> {
    let n_points = rows.len();
    let n_features = rows.first().map(|row| row.len()).unwrap_or(0);

    let mut data = Array2::<f32>::zeros((n_points, n_features));
    for (i, row) in rows.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            data[[i, j]] = *value as f32;
        }
    }

    let config = UmapConfig {
        n_components: settings.n_components,
        manifold: ManifoldParams {
            // flight pins the curve directly, which makes min_dist and spread inert.
            min_dist: 0.0,
            a: Some(settings.curve.a),
            b: Some(settings.curve.b),
            ..Default::default()
        },
        graph: GraphParams {
            n_neighbors: settings.n_neighbours.min(knn.indices.ncols()),
            set_op_mix_ratio: 1.0,
            ..Default::default()
        },
        optimization: OptimizationParams {
            n_epochs: Some(settings.n_epochs),
            ..Default::default()
        },
    };

    let manifold = Umap::new(config.clone()).learn_manifold(data.view(), knn.indices.view(), knn.dists.view());
    let init = spectral_init(manifold.graph(), settings.n_components, settings.seed);

    let mut optimizer = Optimizer::new(
        manifold,
        init,
        settings.n_epochs,
        &config,
        MetricType::Euclidean,
    );
    optimizer.step_epochs(settings.n_epochs, &EuclideanMetric);
    debug!(
        "Embedded {} contigs into {} dimensions with a {} b {}",
        n_points, settings.n_components, settings.curve.a, settings.curve.b
    );
    Ok(optimizer.embedding().mapv(|value| value as f64))
}
