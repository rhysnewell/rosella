use std::collections::HashMap;

use ndarray::{Array2, ArrayView2};
use rand::{Rng, SeedableRng, rngs::StdRng};

use crate::clustering::validity::{dbcv, dbcv_weighted};

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

    fn range(&self) -> ScoreRange;

    fn thresholds(&self) -> ScoreThresholds;
}

/// The interval `score` returns.
#[derive(Debug, Clone, Copy)]
pub struct ScoreRange {
    pub worst: f64,
    pub best: f64,
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

pub const OBJECTIVE_NAMES: [&str; 4] = ["dbcv", "dbcv-bp", "dbcv-floor", "dbcv-bp-floor"];

#[derive(Debug, Clone, Copy, Default)]
pub struct ObjectiveChoice {
    pub weight: ClusterWeight,
    pub floor: bool,
}

impl ObjectiveChoice {
    pub fn parse(name: &str) -> Option<Self> {
        let (weight, floor) = match name {
            "dbcv" => (ClusterWeight::Count, false),
            "dbcv-bp" => (ClusterWeight::Bp, false),
            "dbcv-floor" => (ClusterWeight::Count, true),
            "dbcv-bp-floor" => (ClusterWeight::Bp, true),
            _ => return None,
        };
        Some(Self { weight, floor })
    }

    pub fn build<'a>(&self, lengths: &'a [usize], min_bin_size: usize) -> Dbcv<'a> {
        Dbcv::new(lengths)
            .with_weight(self.weight)
            .with_floor(self.floor.then_some(min_bin_size))
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

    fn range(&self) -> ScoreRange {
        ScoreRange {
            worst: -1.0,
            best: 1.0,
        }
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
