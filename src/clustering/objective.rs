use crate::clustering::codelength::codelength_saving;
use crate::clustering::modularity::modularity;
use crate::embedding::Graph;

/// What a candidate labelling is scored on during the parameter sweep.
pub trait ClusterObjective: Sync {
    fn score_graph(&self, graph: &Graph, labels: &[i32]) -> f64;
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

}
