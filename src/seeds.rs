use rand::{Rng, SeedableRng, rngs::StdRng};

/// One seed per stage that moves a bin. Sampling rides the run seed because holding the rest
/// still and moving it changed nothing on any site.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Seeds {
    pub seed: u64,
    pub knn: u64,
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
