/// One seed per stochastic stage. A single seed across all four makes the spread between
/// runs impossible to attribute, because holding three stages still and moving the fourth
/// is the only thing that separates them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Seeds {
    pub knn: u64,
    pub init: u64,
    pub layout: u64,
    pub sample: u64,
}
