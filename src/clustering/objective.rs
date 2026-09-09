use crate::clustering::codelength::codelength_saving;
use crate::clustering::modularity::modularity;
use crate::embedding::Graph;

/// What a candidate labelling is scored on during the parameter sweep, together with the
/// thresholds that mean nothing away from that score's scale. They travel as one trait
/// because changing the score without re-deriving the thresholds silently breaks refinement.
pub trait ClusterObjective: Sync {
    fn score_graph(&self, graph: &Graph, labels: &[i32]) -> f64;

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

pub const OBJECTIVE_NAMES: [&str; 2] = ["modularity", "codelength"];

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ObjectiveChoice {
    Modularity,
    #[default]
    Codelength,
}

impl ObjectiveChoice {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "modularity" => Some(Self::Modularity),
            "codelength" => Some(Self::Codelength),
            _ => None,
        }
    }

    pub fn build(&self) -> Objective {
        Objective { score: *self }
    }
}

/// Leiden maximises CPM, which carries no degree term, so scoring every rung of its ladder at
/// one fixed modularity resolution ranks them on a criterion the search did not optimise.
pub const EVALUATION_GAMMA: f64 = 1.0;

pub struct Objective {
    score: ObjectiveChoice,
}

impl ClusterObjective for Objective {
    fn score_graph(&self, graph: &Graph, labels: &[i32]) -> f64 {
        match self.score {
            ObjectiveChoice::Modularity => modularity(graph, labels, EVALUATION_GAMMA),
            ObjectiveChoice::Codelength => codelength_saving(graph, labels),
        }
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
