use std::collections::HashMap;

use anyhow::Result;
use rand::{Rng, SeedableRng, rngs::StdRng};

use crate::learned::pieces::Pieces;

pub const CUTS: usize = 4;
pub const SIDES: usize = 2 * CUTS;
const LOW: f64 = 0.3;
const HIGH: f64 = 0.7;

// The two sides of a cut share no k-mer. Overlapping crops let a contig match itself on exact
// composition, and the space then learned nothing about genomes the labels missed.
pub fn pieces(
    assembly: &str,
    wanted: &HashMap<&str, usize>,
    kmer_size: usize,
    seed: u64,
) -> Result<Pieces> {
    crate::learned::pieces::pieces(assembly, wanted, kmer_size, seed, SIDES, sides)
}

pub fn sides(key: u64, length: usize) -> Vec<(usize, usize)> {
    let mut rng = StdRng::seed_from_u64(key);
    (0..CUTS)
        .flat_map(|_| {
            let at = (rng.random_range(LOW..HIGH) * length as f64) as usize;
            [(0, at), (at, length)]
        })
        .collect()
}
