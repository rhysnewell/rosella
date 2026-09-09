use ndarray::Array2;
use umap_rs::{GraphParams, ManifoldParams, Umap, UmapConfig};

use crate::embedding::{
    Graph,
    knn::KnnGraph,
};

const SMALL_DATASET: usize = 10_000;
/// Parameters of UMAP's distance to probability curve, `1 / (1 + a * x^(2b))`.
#[derive(Debug, Clone, Copy)]
pub struct CurveParams {
    pub a: f32,
    pub b: f32,
}

/// `b` was derived from a ratio whose terms are floored so hard it cannot leave a 0.02 wide
/// band, and every value in that band clamps to the ceiling. It is a constant, and no
/// `min_dist` and `spread` fit reaches it either. `a` is not: it reads 1.40 on CAMI I low
/// and 1.47 on high, so its derivation stayed.
const CURVE_B: f32 = 0.6;

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

pub fn curve_params(contig_lengths: &[usize]) -> CurveParams {
    let a = ((n_x(contig_lengths, 10.0).max(1) as f64).log10() * 0.1 + 1.0).clamp(1.4, 2.0);
    CurveParams {
        a: a as f32,
        b: CURVE_B,
    }
}

/// Either flight's pinned pair or the least squares fit UMAP takes from `min_dist` and
/// `spread`, which is the only other way the curve has ever been set.
#[derive(Debug, Clone, Copy)]
pub enum Curve {
    Pinned(CurveParams),
    Fit { min_dist: f32, spread: f32 },
}

impl Curve {
    pub fn from_overrides(contig_lengths: &[usize], overrides: &EmbedOverrides) -> Self {
        match (overrides.min_dist, overrides.spread) {
            (None, None) => {
                let derived = curve_params(contig_lengths);
                Self::Pinned(CurveParams {
                    a: overrides.a.unwrap_or(derived.a),
                    b: overrides.b.unwrap_or(derived.b),
                })
            }
            (min_dist, spread) => Self::Fit {
                min_dist: min_dist.unwrap_or(0.0),
                spread: spread.unwrap_or(1.0),
            },
        }
    }
}

/// Set from the CLI so an ablation can hold the curve or the dimensionality still.
#[derive(Debug, Clone, Copy, Default)]
pub struct EmbedOverrides {
    pub a: Option<f32>,
    pub b: Option<f32>,
    pub min_dist: Option<f32>,
    pub spread: Option<f32>,
    pub knn_candidates: Option<usize>,
    pub graph_weights: crate::embedding::manifold::GraphWeights,
}



pub fn default_epochs(n_points: usize) -> usize {
    if n_points <= SMALL_DATASET { 500 } else { 200 }
}

pub fn manifold_graph(
    n_points: usize,
    knn: &KnnGraph,
    n_neighbours: usize,
    curve: Curve,
) -> (Graph, CurveParams) {
    let config = UmapConfig {
        manifold: match curve {
            Curve::Pinned(curve) => ManifoldParams {
                min_dist: 0.0,
                a: Some(curve.a),
                b: Some(curve.b),
                ..Default::default()
            },
            Curve::Fit { min_dist, spread } => ManifoldParams {
                min_dist,
                spread,
                a: None,
                b: None,
            },
        },
        graph: GraphParams {
            n_neighbors: n_neighbours.min(knn.indices.ncols()),
            set_op_mix_ratio: 1.0,
            ..Default::default()
        },
        ..Default::default()
    };

    let _timer = crate::timing::scope("manifold");
    let empty = Array2::<f32>::zeros((n_points, 0));
    let manifold =
        Umap::new(config).learn_manifold(empty.view(), knn.indices.view(), knn.dists.view());
    let (a, b) = manifold.curve_params();
    (manifold.graph().clone(), CurveParams { a, b })
}

