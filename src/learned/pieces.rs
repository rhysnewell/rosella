use std::collections::HashMap;

use anyhow::{Result, bail};
use ndarray::Array2;
use needletail::Sequence;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::kmers::kmer_counting::{Block, block_of, counts_of, frequencies};

pub const STRETCH_FLOOR: usize = 1000;
pub const VIEWS: usize = 8;

pub struct Pieces {
    pub composition: Vec<f32>,
    pub spans: Vec<f32>,
}

pub fn pieces(
    assembly: &str,
    wanted: &HashMap<&str, usize>,
    kmer_size: usize,
    seed: u64,
) -> Result<Pieces> {
    let block = block_of(kmer_size);
    let width = block.width;
    let mut composition = vec![0.0; wanted.len() * VIEWS * width];
    let mut spans = vec![0.0; wanted.len() * VIEWS];
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
                .map(|(slot, key, sequence)| {
                    Ok((slot, view_rows(&sequence, &block, key, kmer_size)?))
                })
                .collect::<Result<Vec<_>>>()?;
            for (slot, (rows, lengths)) in cut {
                composition[slot * VIEWS * width..(slot + 1) * VIEWS * width]
                    .copy_from_slice(&rows);
                spans[slot * VIEWS..(slot + 1) * VIEWS].copy_from_slice(&lengths);
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

pub fn view_key(seed: u64, name: &str) -> u64 {
    crate::digest::mix(seed ^ crate::digest::fold(name.as_bytes()))
}

// The whole contig, then stretches that start anywhere and run from a floor to the end at most,
// so a long contig is also seen as short as the shortest contig.
pub fn views(key: u64, length: usize) -> Vec<(usize, usize)> {
    let mut rng = StdRng::seed_from_u64(key);
    let floor = STRETCH_FLOOR.min(length);
    std::iter::once((0, length))
        .chain((1..VIEWS).map(|_| {
            let start = rng.random_range(0..=length - floor);
            let end = start + rng.random_range(floor..=length - start);
            (start, end)
        }))
        .collect()
}

fn view_rows(
    sequence: &[u8],
    block: &Block,
    key: u64,
    kmer_size: usize,
) -> Result<(Vec<f32>, Vec<f32>)> {
    let length = sequence.len();
    let mut rows = Vec::with_capacity(VIEWS * block.width);
    let mut lengths = Vec::with_capacity(VIEWS);
    for (start, end) in views(key, length) {
        let (counts, kmers) = counts_of(&sequence[start..end], block);
        rows.extend(frequencies(&counts, kmers));
        lengths.push(end - start);
    }
    let mut table = Array2::from_shape_vec((lengths.len(), block.width), rows)?;
    crate::kmers::clr::clr(&mut table, &lengths, kmer_size)?;
    Ok((
        table.iter().map(|value| *value as f32).collect(),
        lengths.iter().map(|length| *length as f32).collect(),
    ))
}
