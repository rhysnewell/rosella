use ndarray::Array2;

use crate::clustering::validity::dbcv;

/// What a candidate labelling is scored on during the parameter sweep, together with the
/// thresholds that mean nothing away from that score's scale. They travel as one trait
/// because changing the score without re-deriving the thresholds silently breaks refinement.
pub trait ClusterObjective: Sync {
    /// `contigs[i]` is the contig that row `i` of `points` came from. Geometry ignores it;
    /// anything read off the sequence needs it.
    fn score(&self, points: &Array2<f64>, contigs: &[usize], labels: &[i32]) -> f64;

    fn range(&self) -> ScoreRange;

    fn thresholds(&self) -> ScoreThresholds;
}

/// The interval `score` returns.
#[derive(Debug, Clone, Copy)]
pub struct ScoreRange {
    pub worst: f64,
    pub best: f64,
}

/// The four refinement decisions read off the objective's scale rather than off the data.
#[derive(Debug, Clone, Copy)]
pub struct ScoreThresholds {
    /// Above this the first clustering is good enough that re-embedding cannot improve it.
    pub re_embed_ceiling: f64,
    /// A split into fewer than two clusters has to reach this to be believed.
    pub single_cluster: f64,
    /// The bar for a bin that tripped one of the distance levels.
    pub tripped_ceiling: f64,
    /// The bar for a bin that is merely dirty.
    pub dirty_ceiling: f64,
}

/// Density based cluster validity. Asks how well separated a labelling is, which is not the
/// same question as how much each cluster looks like a genome.
pub struct Dbcv;

impl ClusterObjective for Dbcv {
    fn score(&self, points: &Array2<f64>, _contigs: &[usize], labels: &[i32]) -> f64 {
        dbcv(points, labels)
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
            tripped_ceiling: 0.5,
            dirty_ceiling: 0.95,
        }
    }
}
