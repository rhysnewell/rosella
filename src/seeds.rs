use rand::{Rng, SeedableRng, rngs::StdRng};

/// One seed per stochastic stage. A single seed across all of them makes the spread between
/// runs impossible to attribute, because holding the rest still and moving one is the only
/// thing that separates them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Seeds {
    pub knn: u64,
    pub sample: u64,
    pub partition: u64,
}

/// Partial Fisher-Yates: the first `k` swaps of a full shuffle, so drawing a sample costs the
/// sample and not the population.
pub fn sample_positions(n: usize, k: usize, seed: u64) -> Vec<usize> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut positions = (0..n).collect::<Vec<_>>();
    for position in 0..k {
        positions.swap(position, rng.random_range(position..n));
    }
    positions.truncate(k);
    positions
}
