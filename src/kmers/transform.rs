use anyhow::{Result, bail};
use ndarray::Array2;
use rayon::prelude::*;

use crate::kmers::kmer_counting::canonical_index;

const BASES: &[u8; 4] = b"ACGT";

pub fn sqrt_frequencies(table: &mut Array2<f64>) {
    let n_cols = table.ncols();
    for mut row in table.rows_mut() {
        let sum = row.iter().filter(|value| **value > 0.0).sum::<f64>();
        let scale = if sum > 0.0 {
            1.0 / sum
        } else {
            1.0 / n_cols as f64
        };
        row.map_inplace(|value| *value = (value.max(0.0) * scale).sqrt());
    }
}

/// Observed against a maximal-order Markov expectation, Teeling et al. 2004. The z score is
/// constant on a canonical pair, so the table keeps its columns.
pub fn tetra_z(table: &mut Array2<f64>, kmer_size: usize) -> Result<()> {
    if kmer_size < 3 {
        bail!("tetra-z needs a k-mer of at least 3 bases, got {kmer_size}.");
    }
    let columns = column_of_code(kmer_size, table.ncols())?;
    let representatives = representatives(&columns, table.ncols(), kmer_size);

    let rows = table
        .rows()
        .into_iter()
        .collect::<Vec<_>>()
        .into_par_iter()
        .flat_map(|row| {
            let full = columns
                .iter()
                .map(|column| row[*column])
                .collect::<Vec<_>>();
            let (prefix, suffix, middle) = marginals(&full, kmer_size);
            representatives
                .iter()
                .map(|(code, prefix_code, suffix_code, middle_code)| {
                    z_score(
                        full[*code],
                        prefix[*prefix_code],
                        suffix[*suffix_code],
                        middle[*middle_code],
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    *table = Array2::from_shape_vec(table.dim(), rows)?;
    Ok(())
}

fn z_score(observed: f64, prefix: f64, suffix: f64, middle: f64) -> f64 {
    if middle <= 0.0 {
        return 0.0;
    }
    let expected = prefix * suffix / middle;
    let variance = expected * (middle - prefix) * (middle - suffix) / (middle * middle);
    if variance <= 0.0 {
        return 0.0;
    }
    (observed - expected) / variance.sqrt()
}

fn marginals(full: &[f64], kmer_size: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let inner = 4usize.pow(kmer_size as u32 - 1);
    let core = 4usize.pow(kmer_size as u32 - 2);
    let mut prefix = vec![0.0; inner];
    let mut suffix = vec![0.0; inner];
    let mut middle = vec![0.0; core];
    for (code, value) in full.iter().enumerate() {
        prefix[code / 4] += value;
        suffix[code % inner] += value;
        middle[(code / 4) % core] += value;
    }
    (prefix, suffix, middle)
}

/// Column of the canonical pair each of the `4^k` k-mers belongs to, indexed by base-4 code.
fn column_of_code(kmer_size: usize, n_cols: usize) -> Result<Vec<usize>> {
    let index = canonical_index(kmer_size);
    if index.len() != n_cols {
        bail!(
            "A k-mer table of {n_cols} columns does not match the {} canonical {kmer_size}-mers.",
            index.len()
        );
    }
    (0..4usize.pow(kmer_size as u32))
        .map(|code| {
            let kmer = kmer_of(code, kmer_size);
            let reverse = needletail::Sequence::reverse_complement(&kmer[..]);
            let canonical = if kmer <= reverse { &kmer } else { &reverse };
            index
                .get(canonical)
                .copied()
                .ok_or_else(|| anyhow::anyhow!("Canonical k-mer missing from the table index."))
        })
        .collect()
}

/// One k-mer code per column, with the sub-k-mer codes its expectation reads.
fn representatives(
    columns: &[usize],
    n_cols: usize,
    kmer_size: usize,
) -> Vec<(usize, usize, usize, usize)> {
    let inner = 4usize.pow(kmer_size as u32 - 1);
    let core = 4usize.pow(kmer_size as u32 - 2);
    let mut chosen = vec![None; n_cols];
    for (code, column) in columns.iter().enumerate() {
        chosen[*column].get_or_insert(code);
    }
    chosen
        .into_iter()
        .map(|code| {
            let code = code.expect("every column has at least one k-mer");
            (code, code / 4, code % inner, (code / 4) % core)
        })
        .collect()
}

fn kmer_of(code: usize, kmer_size: usize) -> Vec<u8> {
    (0..kmer_size)
        .map(|position| {
            let shift = 2 * (kmer_size - 1 - position);
            BASES[(code >> shift) & 3]
        })
        .collect()
}
