use anyhow::Result;
use log::debug;
use ndarray::Array2;
use umap_rs::{GraphParams, ManifoldParams, OptimizationParams, Umap, UmapConfig};

use crate::embedding::{
    knn::KnnGraph,
    layout::{LayoutSettings, optimise},
    spectral::spectral_init,
};
use crate::seeds::Seeds;

const MIN_COMPONENTS: usize = 2;
const SMALL_DATASET: usize = 10_000;
const MAX_COMPONENTS: usize = 10;

/// How far a single contig's length is allowed to move its edge sampling rate. Unbounded,
/// a megabase contig would be drawn a thousand times more often than a 1.5 kb one and the
/// short contigs would never move.
const LENGTH_WEIGHT_RANGE: (f32, f32) = (0.25, 4.0);

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
///
/// The floors on both terms hold the ratio inside [0.82, 1.0], so `b` is a constant wearing
/// a derivation and the window is what decides it. CAMI I low and medium both rank 0.6 above
/// 0.5 above 0.4 above flight's 0.3, and low collapses at 0.7, so the window ends at 0.6.
pub fn curve_params(contig_lengths: &[usize]) -> CurveParams {
    let log10 = |value: usize| (value.max(1) as f64).log10();

    let numerator = log10(n_x(contig_lengths, 25.0)).max(50_000f64.log10());
    let denominator = 500_000f64.log10().max(log10(n_x(contig_lengths, 75.0)));
    let b = (0.1 * (numerator / denominator) + 0.52).clamp(0.5, 0.6);
    let a = (log10(n_x(contig_lengths, 10.0)) * 0.1 + 1.0).clamp(1.4, 2.0);

    CurveParams {
        a: a as f32,
        b: b as f32,
    }
}

/// The sample count was the wrong input: it gave 2 dimensions to a single-sample assembly
/// whose data occupies 7, and a graph with more near-disconnected groups than dimensions
/// leaves the spectral start undetermined. Falls back to the old rule when the estimate
/// cannot be taken.
pub fn n_components(intrinsic_dimension: Option<f64>, n_samples: usize) -> usize {
    let dimensions = match intrinsic_dimension {
        Some(estimate) if estimate.is_finite() && estimate >= 1.0 => estimate.round() as usize,
        _ => n_samples,
    };
    dimensions.clamp(MIN_COMPONENTS, MAX_COMPONENTS)
}

/// Set from the CLI to pin what the assembly would otherwise decide, so an ablation can
/// hold the curve or the dimensionality still.
#[derive(Debug, Clone, Copy, Default)]
pub struct EmbedOverrides {
    pub a: Option<f32>,
    pub b: Option<f32>,
    pub n_components: Option<usize>,
    pub n_epochs: Option<usize>,
    pub length_weight: f64,
}

pub struct EmbedSettings {
    pub n_components: usize,
    pub n_neighbours: usize,
    pub curve: CurveParams,
    pub n_epochs: usize,
    pub seeds: Seeds,
    pub vertex_weights: Vec<f32>,
}

/// Per-contig edge sampling weights, empty at power 0 so the layout is untouched. UMAP has
/// no per-point weight, so length enters through how often a contig's edges are drawn, and
/// the geometric mean of 1 leaves the threshold that drops the weakest edges where it was.
pub fn length_weights(lengths: &[usize], power: f64) -> Vec<f32> {
    if power == 0.0 || lengths.is_empty() {
        return Vec::new();
    }

    let mut sorted = lengths.to_vec();
    sorted.sort_unstable();
    let reference = sorted[sorted.len() / 2].max(1) as f64;

    let weights = lengths
        .iter()
        .map(|length| {
            ((*length as f64 / reference).powf(power) as f32)
                .clamp(LENGTH_WEIGHT_RANGE.0, LENGTH_WEIGHT_RANGE.1)
        })
        .collect::<Vec<f32>>();

    let log_mean = weights.iter().map(|w| (*w as f64).ln()).sum::<f64>() / weights.len() as f64;
    let scale = log_mean.exp() as f32;
    weights.into_iter().map(|w| w / scale).collect()
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

    let manifold = {
        let _timer = crate::timing::scope("manifold");
        Umap::new(config).learn_manifold(data.view(), knn.indices.view(), knn.dists.view())
    };
    let init = {
        let _timer = crate::timing::scope("spectral_init");
        spectral_init(manifold.graph(), settings.n_components, settings.seeds.init)
    };

    let layout = LayoutSettings {
        curve: settings.curve,
        n_epochs: settings.n_epochs,
        seed: settings.seeds.layout,
    };
    let embedding = {
        let _timer = crate::timing::scope("layout_sgd");
        optimise(manifold.graph(), init, &layout, &settings.vertex_weights)
    };
    debug!(
        "Embedded {} contigs into {} dimensions with a {} b {}",
        n_points, settings.n_components, settings.curve.a, settings.curve.b
    );
    Ok(embedding.mapv(|value| value as f64))
}
