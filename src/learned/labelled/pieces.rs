use std::collections::HashMap;

use anyhow::Result;
use rand::{Rng, SeedableRng, rngs::StdRng};

use crate::learned::pieces::{Pieces, STRETCH_FLOOR};

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

fn sides(key: u64, length: usize) -> Vec<(usize, usize)> {
    let mut rng = StdRng::seed_from_u64(key);
    (0..CUTS)
        .flat_map(|_| {
            let at = (rng.random_range(LOW..HIGH) * length as f64) as usize;
            [stretch(&mut rng, 0, at), stretch(&mut rng, at, length)]
        })
        .collect()
}

// A half never shows the composition noise of a short stretch, so a genome that varies along
// its length, as a eukaryote does, reads as several.
fn stretch(rng: &mut StdRng, start: usize, end: usize) -> (usize, usize) {
    let size = end - start;
    let span = rng.random_range(STRETCH_FLOOR.min(size)..=size);
    let from = start + rng.random_range(0..=size - span);
    (from, from + span)
}
