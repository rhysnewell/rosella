use std::collections::HashMap;

use ndarray::{Array2, ArrayView2};
use rand::{Rng, SeedableRng, rngs::StdRng};

use crate::clustering::codelength::codelength_saving;
use crate::clustering::modularity::modularity;
use crate::clustering::validity::{dbcv, dbcv_weighted};
use crate::embedding::Graph;

/// Validity is quadratic in the points it scores and the sweep scores every combination, so
/// a large embedding is scored on a sample of itself.
const VALIDITY_SAMPLE_LIMIT: usize = 5000;

/// What a candidate labelling is scored on during the parameter sweep, together with the
/// thresholds that mean nothing away from that score's scale. They travel as one trait
/// because changing the score without re-deriving the thresholds silently breaks refinement.
pub trait ClusterObjective: Sync {
    /// `labels` and `contigs` cover every row while `sample` is the subset the geometry is
    /// measured on, because a whole cluster keeps only a fraction of its contigs in a sample.
    fn score(&self, sample: &EmbeddingSample, contigs: &[usize], labels: &[i32]) -> f64;

    /// `None` means this objective has nothing to say about a graph, so the caller falls
    /// back to `score`.
    fn score_graph(&self, _graph: &Graph, _labels: &[i32]) -> Option<f64> {
        None
    }

    /// Asked before the embedding is built rather than after, so a caller that would only
    /// have fed the layout to `score_graph` can skip building one.
    fn needs_layout(&self) -> bool {
        true
    }

    fn thresholds(&self) -> ScoreThresholds;
}

/// The refinement decisions read off the objective's scale. The bar a split has to clear is
/// not among them: it is in distance units and lives in `refine::bar`.
#[derive(Debug, Clone, Copy)]
pub struct ScoreThresholds {
    /// Above this the first clustering is good enough that re-embedding cannot improve it.
    pub re_embed_ceiling: f64,
    /// A split into fewer than two clusters has to reach this to be believed.
    pub single_cluster: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ClusterWeight {
    /// Contigs, which is what flight counts.
    #[default]
    Count,
    /// Assembled bases, which is what the bin floor is applied in.
    Bp,
}

pub const OBJECTIVE_NAMES: [&str; 6] = [
    "dbcv",
    "dbcv-bp",
    "dbcv-floor",
    "dbcv-bp-floor",
    "modularity",
    "codelength",
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GraphScore {
    Modularity,
    Codelength,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct ObjectiveChoice {
    pub weight: ClusterWeight,
    pub floor: bool,
    pub graph: Option<GraphScore>,
}

impl ObjectiveChoice {
    pub fn parse(name: &str) -> Option<Self> {
        let (weight, floor, graph) = match name {
            "dbcv" => (ClusterWeight::Count, false, None),
            "dbcv-bp" => (ClusterWeight::Bp, false, None),
            "dbcv-floor" => (ClusterWeight::Count, true, None),
            "dbcv-bp-floor" => (ClusterWeight::Bp, true, None),
            "modularity" => (ClusterWeight::Count, false, Some(GraphScore::Modularity)),
            "codelength" => (ClusterWeight::Count, false, Some(GraphScore::Codelength)),
            _ => return None,
        };
        Some(Self {
            weight,
            floor,
            graph,
        })
    }

    pub fn build<'a>(&self, lengths: &'a [usize], min_bin_size: usize) -> Objective<'a> {
        let layout = Dbcv::new(lengths)
            .with_weight(self.weight)
            .with_floor(self.floor.then_some(min_bin_size));
        match self.graph {
            Some(score) => Objective::Graph(GraphObjective { layout, score }),
            None => Objective::Dbcv(layout),
        }
    }
}

pub enum Objective<'a> {
    Dbcv(Dbcv<'a>),
    Graph(GraphObjective<'a>),
}

impl ClusterObjective for Objective<'_> {
    fn score(&self, sample: &EmbeddingSample, contigs: &[usize], labels: &[i32]) -> f64 {
        match self {
            Self::Dbcv(inner) => inner.score(sample, contigs, labels),
            Self::Graph(inner) => inner.score(sample, contigs, labels),
        }
    }

    fn score_graph(&self, graph: &Graph, labels: &[i32]) -> Option<f64> {
        match self {
            Self::Dbcv(inner) => inner.score_graph(graph, labels),
            Self::Graph(inner) => inner.score_graph(graph, labels),
        }
    }

    fn needs_layout(&self) -> bool {
        match self {
            Self::Dbcv(inner) => inner.needs_layout(),
            Self::Graph(inner) => inner.needs_layout(),
        }
    }

    fn thresholds(&self) -> ScoreThresholds {
        match self {
            Self::Dbcv(inner) => inner.thresholds(),
            Self::Graph(inner) => inner.thresholds(),
        }
    }
}

/// Leiden maximises CPM, which carries no degree term, so scoring every rung of its ladder at
/// one fixed modularity resolution ranks them on a criterion the search did not optimise.
pub const EVALUATION_GAMMA: f64 = 1.0;

/// Falls back to DBCV wherever no graph exists, which is the HDBSCAN sweep and the
/// refiner's in-place attempt.
pub struct GraphObjective<'a> {
    layout: Dbcv<'a>,
    score: GraphScore,
}

impl ClusterObjective for GraphObjective<'_> {
    fn score(&self, sample: &EmbeddingSample, contigs: &[usize], labels: &[i32]) -> f64 {
        self.layout.score(sample, contigs, labels)
    }

    fn score_graph(&self, graph: &Graph, labels: &[i32]) -> Option<f64> {
        Some(match self.score {
            GraphScore::Modularity => modularity(graph, labels, EVALUATION_GAMMA),
            GraphScore::Codelength => codelength_saving(graph, labels),
        })
    }

    fn needs_layout(&self) -> bool {
        false
    }

    /// Both scores put one community at exactly zero and rank descending, so they share these
    /// to keep the ladder ranking the only difference between them.
    fn thresholds(&self) -> ScoreThresholds {
        ScoreThresholds {
            re_embed_ceiling: 0.85,
            single_cluster: 0.7,
        }
    }
}

/// Asks how well separated a labelling is, not how much each cluster looks like a genome.
/// `weight` and `floor` exist because counting contigs ranks in a unit the writer never
/// uses, and both are off by default.
pub struct Dbcv<'a> {
    weight: ClusterWeight,
    floor: Option<usize>,
    lengths: &'a [usize],
}

impl<'a> Dbcv<'a> {
    pub fn new(lengths: &'a [usize]) -> Self {
        Self {
            weight: ClusterWeight::Count,
            floor: None,
            lengths,
        }
    }

    pub fn with_weight(mut self, weight: ClusterWeight) -> Self {
        self.weight = weight;
        self
    }

    /// Clusters under `floor` assembled bases score nothing, because the writer is going to
    /// dissolve them whatever the geometry says.
    pub fn with_floor(mut self, floor: Option<usize>) -> Self {
        self.floor = floor;
        self
    }

    fn is_default(&self) -> bool {
        self.weight == ClusterWeight::Count && self.floor.is_none()
    }

    fn bp_of(&self, contig: usize) -> usize {
        self.lengths.get(contig).copied().unwrap_or(0)
    }
}

impl ClusterObjective for Dbcv<'_> {
    fn score(&self, sample: &EmbeddingSample, contigs: &[usize], labels: &[i32]) -> f64 {
        let sampled_labels = sample.labels(labels);

        if self.is_default() {
            return dbcv(&sample.rows, &sampled_labels);
        }

        let mut bp: HashMap<i32, usize> = HashMap::new();
        let mut total = 0usize;
        for (index, label) in labels.iter().enumerate() {
            let length = self.bp_of(contigs[index]);
            total += length;
            if *label >= 0 {
                *bp.entry(*label).or_default() += length;
            }
        }

        let total = total as f64;
        let scored = sampled_labels.len() as f64;
        dbcv_weighted(&sample.rows, &sampled_labels, |label, size| {
            let cluster_bp = bp.get(&label).copied().unwrap_or(0);
            if self.floor.is_some_and(|floor| cluster_bp < floor) {
                return 0.0;
            }
            match self.weight {
                ClusterWeight::Count => size as f64 / scored,
                ClusterWeight::Bp if total > 0.0 => cluster_bp as f64 / total,
                ClusterWeight::Bp => 0.0,
            }
        })
    }

    fn thresholds(&self) -> ScoreThresholds {
        ScoreThresholds {
            re_embed_ceiling: 0.95,
            single_cluster: 0.9,
        }
    }
}

/// Built once per sweep rather than per parameter combination: the sweep scores 55 of them
/// against one embedding, and rebuilding this for each cost four times the clustering stage.
pub struct EmbeddingSample {
    indices: Vec<usize>,
    rows: Array2<f64>,
}

impl EmbeddingSample {
    pub fn new(points: ArrayView2<f64>, seed: u64) -> Self {
        let indices = validity_sample(points.nrows(), seed);
        let rows = sample_rows(&points, &indices);
        Self { indices, rows }
    }

    fn labels(&self, all: &[i32]) -> Vec<i32> {
        self.indices.iter().map(|index| all[*index]).collect()
    }
}

fn validity_sample(n: usize, seed: u64) -> Vec<usize> {
    if n <= VALIDITY_SAMPLE_LIMIT {
        return (0..n).collect();
    }

    let mut rng = StdRng::seed_from_u64(seed);
    let mut chosen = (0..n).collect::<Vec<_>>();
    for position in 0..VALIDITY_SAMPLE_LIMIT {
        chosen.swap(position, rng.random_range(position..n));
    }
    chosen.truncate(VALIDITY_SAMPLE_LIMIT);
    chosen.sort_unstable();
    chosen
}

fn sample_rows(points: &ArrayView2<f64>, sample: &[usize]) -> Array2<f64> {
    let mut rows = Array2::zeros((sample.len(), points.ncols()));
    for (position, index) in sample.iter().enumerate() {
        rows.row_mut(position).assign(&points.row(*index));
    }
    rows
}
