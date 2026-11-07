use rand::{SeedableRng, rngs::StdRng};

// One seed per stage that moves a bin. Sampling rides the run seed because holding the rest
// still and moving it changed nothing on any site.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Seeds {
    pub seed: u64,
    pub knn: u64,
    pub partition: u64,
}

pub fn sample_positions(n: usize, k: usize, seed: u64) -> Vec<usize> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut positions = (0..n).collect::<Vec<_>>();
    let k = k.min(n);
    for at in 0..k {
        crate::rows::draw(&mut rng, &mut positions, at);
    }
    positions.truncate(k);
    positions
}
