use std::collections::HashMap;

use anyhow::{Result, bail};
use ndarray::Array2;
use needletail::Sequence;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::kmers::kmer_counting::{Block, block_of, counts_of, frequencies};
use crate::learned::pieces::{Pieces, view_key};

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
    let block = block_of(kmer_size);
    let width = block.width;
    let mut composition = vec![0.0; wanted.len() * SIDES * width];
    let mut spans = vec![0.0; wanted.len() * SIDES];
    let mut found = 0;
    crate::kmers::pipelined(
        assembly,
        |record| {
            let name = crate::contig_id(record.id())?;
            Ok(wanted.get(name).map(|slot| {
                (
                    *slot,
                    view_key(seed, name),
                    record.normalize(false).into_owned(),
                )
            }))
        },
        |chunk| {
            let cut = chunk
                .into_par_iter()
                .flatten()
                .map(|(slot, key, sequence)| Ok((slot, sides(&sequence, &block, key, kmer_size)?)))
                .collect::<Result<Vec<_>>>()?;
            for (slot, (rows, lengths)) in cut {
                composition[slot * SIDES * width..(slot + 1) * SIDES * width]
                    .copy_from_slice(&rows);
                spans[slot * SIDES..(slot + 1) * SIDES].copy_from_slice(&lengths);
                found += 1;
            }
            Ok(())
        },
    )?;
    if found != wanted.len() {
        bail!(
            "{} of {} contigs to cut were found in {assembly}",
            found,
            wanted.len()
        );
    }
    Ok(Pieces { composition, spans })
}

fn sides(
    sequence: &[u8],
    block: &Block,
    key: u64,
    kmer_size: usize,
) -> Result<(Vec<f32>, Vec<f32>)> {
    let mut rng = StdRng::seed_from_u64(key);
    let length = sequence.len();
    let mut rows = Vec::with_capacity(SIDES * block.width);
    let mut lengths = Vec::with_capacity(SIDES);
    for _ in 0..CUTS {
        let at = (rng.random_range(LOW..HIGH) * length as f64) as usize;
        for side in [&sequence[..at], &sequence[at..]] {
            let (counts, kmers) = counts_of(side, block);
            rows.extend(frequencies(&counts, kmers));
            lengths.push(side.len());
        }
    }
    let mut table = Array2::from_shape_vec((SIDES, block.width), rows)?;
    crate::kmers::clr::clr(&mut table, &lengths, kmer_size)?;
    Ok((
        table.iter().map(|value| *value as f32).collect(),
        lengths.iter().map(|length| *length as f32).collect(),
    ))
}
