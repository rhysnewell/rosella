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

pub fn mix(value: u64) -> u64 {
    let mut mixed = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    mixed = (mixed ^ (mixed >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    mixed ^ (mixed >> 31)
}

// Each item keeps one key per seed, so two pools that share most items share most of the sample
// rather than drawing it again. Items come back in their input order.
pub fn consistent_sample(items: &[usize], k: usize, seed: u64) -> Vec<usize> {
    if items.len() <= k {
        return items.to_vec();
    }
    let key = |item: usize| mix(mix(item as u64) ^ seed);
    let mut keys = items.iter().map(|item| key(*item)).collect::<Vec<_>>();
    let bar = *keys.select_nth_unstable(k - 1).1;
    items
        .iter()
        .copied()
        .filter(|item| key(*item) <= bar)
        .collect()
}
