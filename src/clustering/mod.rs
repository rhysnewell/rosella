pub mod cluster_utils;
pub mod clusterer;

pub struct PropagatedLabel {
    pub label: Option<usize>,
    pub weight: f64,
    pub index: usize
}