/// One seed per stochastic stage. A single seed across all of them makes the spread between
/// runs impossible to attribute, because holding the rest still and moving one is the only
/// thing that separates them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Seeds {
    pub knn: u64,
    pub sample: u64,
    pub partition: u64,
}
