use std::collections::{BTreeMap, HashMap};

use rayon::prelude::*;

use crate::kmers::sketch::ContigSketches;
use crate::markers::replicon::LEAST_BASES;

// A genome's pieces share no sequence, so a contig its bin-mates mostly hold again is one more
// copy of an element the assembly kept several times, not a piece of the genome.
const COPY_SHARE: f64 = 0.5;

// Trial switch for the floor arm, removed once one floor ships.
pub fn copy_floor() -> usize {
    std::env::var("ROSELLA_COPY_FLOOR")
        .ok()
        .and_then(|floor| floor.parse().ok())
        .unwrap_or(LEAST_BASES)
}

pub fn copies(
    bins: &BTreeMap<usize, Vec<usize>>,
    sketches: &ContigSketches,
    lengths: &[usize],
    least: usize,
) -> Vec<usize> {
    let mut found = bins
        .par_iter()
        .flat_map_iter(|(_, members)| copied(members, sketches, lengths, least))
        .collect::<Vec<_>>();
    found.sort_unstable();
    found
}

fn copied(
    members: &[usize],
    sketches: &ContigSketches,
    lengths: &[usize],
    least: usize,
) -> Vec<usize> {
    let long = members
        .iter()
        .copied()
        .filter(|contig| lengths[*contig] >= least)
        .collect::<Vec<_>>();
    let mut holders = HashMap::<u64, u32>::new();
    for contig in &long {
        for hash in sketches.sketch(*contig) {
            *holders.entry(*hash).or_default() += 1;
        }
    }
    long.into_iter()
        .filter(|contig| {
            let hashes = sketches.sketch(*contig);
            let shared = hashes.iter().filter(|hash| holders[*hash] > 1).count();
            !hashes.is_empty() && shared as f64 > COPY_SHARE * hashes.len() as f64
        })
        .collect()
}
