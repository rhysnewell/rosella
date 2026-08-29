use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};

use crate::embedding::{Graph, umap::CurveParams};

const GAMMA: f32 = 1.0;
const INITIAL_ALPHA: f32 = 1.0;
const NEGATIVE_SAMPLE_RATE: usize = 5;
const CLIP: f32 = 4.0;

pub struct LayoutSettings {
    pub curve: CurveParams,
    pub n_epochs: usize,
    pub seed: u64,
}

/// UMAP's layout optimisation, replacing umap-rs's. Its version draws negative samples
/// from an OS seeded thread RNG and writes the embedding from many threads without
/// synchronisation, so two runs of identical code at one seed give different bins.
/// This one walks the edges in index order from one seeded stream, which costs the
/// parallelism and buys a reproducible result.
pub fn optimise(
    graph: &Graph,
    mut embedding: Array2<f32>,
    settings: &LayoutSettings,
    vertex_weights: &[f32],
) -> Array2<f32> {
    let edges = Edges::from_graph(graph, settings.n_epochs, vertex_weights);
    if edges.head.is_empty() {
        return embedding;
    }

    let n_vertices = embedding.nrows();
    let dim = embedding.ncols();
    let a = settings.curve.a;
    let b = settings.curve.b;

    let mut next_sample = edges.epochs_per_sample.clone();
    let epochs_per_negative: Vec<f64> = edges
        .epochs_per_sample
        .iter()
        .map(|value| value / NEGATIVE_SAMPLE_RATE as f64)
        .collect();
    let mut next_negative = epochs_per_negative.clone();

    let mut rng = StdRng::seed_from_u64(settings.seed);

    for epoch in 0..settings.n_epochs {
        let alpha = INITIAL_ALPHA * (1.0 - epoch as f32 / settings.n_epochs as f32);
        let now = epoch as f64;

        for edge in 0..edges.head.len() {
            if next_sample[edge] > now {
                continue;
            }

            let head = edges.head[edge] as usize;
            let tail = edges.tail[edge] as usize;
            let grad = attractive_gradient(&embedding, head, tail, dim, a, b);
            for d in 0..dim {
                let step = clip(grad * (embedding[[head, d]] - embedding[[tail, d]])) * alpha;
                embedding[[head, d]] += step;
                embedding[[tail, d]] -= step;
            }
            next_sample[edge] += edges.epochs_per_sample[edge];

            let n_negative = ((now - next_negative[edge]) / epochs_per_negative[edge]) as usize;
            for _ in 0..n_negative {
                let other = rng.random_range(0..n_vertices);
                if other == head {
                    continue;
                }
                let grad = repulsive_gradient(&embedding, head, other, dim, a, b);
                for d in 0..dim {
                    let difference = embedding[[head, d]] - embedding[[other, d]];
                    let step = if grad > 0.0 {
                        clip(grad * difference)
                    } else {
                        CLIP
                    };
                    embedding[[head, d]] += step * alpha;
                }
            }
            next_negative[edge] += n_negative as f64 * epochs_per_negative[edge];
        }
    }

    embedding
}

/// The 1-simplices worth sampling, in row order. Edges too weak to be drawn even once
/// over the whole run are dropped, which is what umap-learn's threshold does.
struct Edges {
    head: Vec<u32>,
    tail: Vec<u32>,
    epochs_per_sample: Vec<f64>,
}

impl Edges {
    /// `vertex_weights` scales how often each edge is drawn, empty meaning uniform. The
    /// maximum is taken over the scaled weights so the drop threshold moves with them,
    /// rather than pruning whatever the scaling pushed down.
    fn from_graph(graph: &Graph, n_epochs: usize, vertex_weights: &[f32]) -> Self {
        let scale = |row: usize, column: u32, weight: f32| {
            if vertex_weights.is_empty() {
                weight
            } else {
                weight * (vertex_weights[row] * vertex_weights[column as usize]).sqrt()
            }
        };

        let mut max_weight = 0.0f32;
        for row in 0..graph.rows() {
            let start = graph.indptr().index(row);
            let end = graph.indptr().index(row + 1);
            for entry in start..end {
                max_weight =
                    max_weight.max(scale(row, graph.indices()[entry], graph.data()[entry]));
            }
        }
        let threshold = max_weight / n_epochs.max(1) as f32;
        let max_weight = if max_weight <= 0.0 { 1.0 } else { max_weight };

        let mut head = Vec::new();
        let mut tail = Vec::new();
        let mut epochs_per_sample = Vec::new();

        for row in 0..graph.rows() {
            let start = graph.indptr().index(row);
            let end = graph.indptr().index(row + 1);
            for entry in start..end {
                let column = graph.indices()[entry];
                let weight = scale(row, column, graph.data()[entry]);
                if weight < threshold {
                    continue;
                }
                head.push(row as u32);
                tail.push(column);
                epochs_per_sample.push((max_weight / weight) as f64);
            }
        }

        Self {
            head,
            tail,
            epochs_per_sample,
        }
    }
}

fn squared_distance(embedding: &Array2<f32>, i: usize, j: usize, dim: usize) -> f32 {
    let mut total = 0.0;
    for d in 0..dim {
        let difference = embedding[[i, d]] - embedding[[j, d]];
        total += difference * difference;
    }
    total
}

fn attractive_gradient(
    embedding: &Array2<f32>,
    head: usize,
    tail: usize,
    dim: usize,
    a: f32,
    b: f32,
) -> f32 {
    let distance = squared_distance(embedding, head, tail, dim);
    if distance <= 0.0 {
        return 0.0;
    }
    let powered = distance.powf(b);
    -2.0 * a * b * powered / distance / (a * powered * distance + 1.0)
}

fn repulsive_gradient(
    embedding: &Array2<f32>,
    head: usize,
    other: usize,
    dim: usize,
    a: f32,
    b: f32,
) -> f32 {
    let distance = squared_distance(embedding, head, other, dim);
    if distance <= 0.0 {
        return 0.0;
    }
    2.0 * GAMMA * b / ((0.001 + distance) * (a * distance.powf(b) + 1.0))
}

fn clip(value: f32) -> f32 {
    value.clamp(-CLIP, CLIP)
}
