use anyhow::{Result, bail};
use log::debug;
use ndarray::Array2;
use rayon::prelude::*;

use crate::kmers::kmer_counting::canonical_count;

// Martin-Fernandez et al. (2003) take 0.65 of the detection limit.
const REPLACEMENT_FRACTION: f64 = 0.65;

pub fn detection_limit(contig_length: usize, kmer_size: usize) -> f64 {
    let positions = contig_length.saturating_sub(kmer_size - 1).max(1);
    REPLACEMENT_FRACTION / positions as f64
}

fn median_replacement(contig_lengths: &[usize], kmer_size: usize) -> f64 {
    if contig_lengths.is_empty() {
        return f64::NAN;
    }
    let mut deltas = contig_lengths
        .iter()
        .map(|length| detection_limit(*length, kmer_size))
        .collect::<Vec<_>>();
    deltas.sort_by(f64::total_cmp);
    deltas[deltas.len() / 2]
}

pub fn clr(table: &mut Array2<f64>, contig_lengths: &[usize], kmer_size: usize) -> Result<()> {
    let n_rows = table.nrows();
    let n_cols = table.ncols();
    if contig_lengths.len() != n_rows {
        bail!(
            "Centre log ratio needs one length per row, got {} lengths for {} contigs.",
            contig_lengths.len(),
            n_rows
        );
    }

    let width = canonical_count(kmer_size);
    if width != n_cols {
        bail!("k-mer size {kmer_size} wants {width} columns and the table holds {n_cols}.");
    }

    debug!(
        "Composition zeros: {:.4} of {} cells (median replacement {:.3e}, k={})",
        table.iter().filter(|value| **value <= 0.0).count() as f64 / (n_rows * n_cols) as f64,
        n_rows * n_cols,
        median_replacement(contig_lengths, kmer_size),
        kmer_size
    );

    if !table.is_standard_layout() {
        *table = table.as_standard_layout().into_owned();
    }
    table
        .as_slice_mut()
        .expect("a standard layout table is contiguous")
        .par_chunks_mut(n_cols.max(1))
        .zip(contig_lengths)
        .for_each(|(row, length)| centre(row, detection_limit(*length, kmer_size)));
    Ok(())
}

fn centre(row: &mut [f64], delta: f64) {
    // ndarray folds eight lanes at once, and a sequential sum would move the table by an ulp.
    let row_sum = ndarray::ArrayView1::from(&*row).sum();
    let n_zeros = row.iter().filter(|value| **value <= 0.0).count();
    let retained = (1.0 - n_zeros as f64 * delta).max(f64::MIN_POSITIVE);
    for value in row.iter_mut() {
        let replaced = if *value <= 0.0 || row_sum <= 0.0 {
            delta
        } else {
            *value / row_sum * retained
        };
        *value = replaced.ln();
    }
    let log_geometric_mean = row.iter().sum::<f64>() / row.len() as f64;
    for value in row.iter_mut() {
        *value -= log_geometric_mean;
    }
}
