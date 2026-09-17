use anyhow::{Result, bail};
use log::debug;
use ndarray::{Array, Array2};
use rayon::prelude::*;

use crate::kmers::kmer_counting::canonical_count;

/// Martin-Fernandez et al. (2003) take 0.65 of the detection limit.
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

/// Each k is its own composition, so the replacement and the geometric mean run per block.
pub fn clr(table: &mut Array2<f64>, contig_lengths: &[usize], kmer_sizes: &[usize]) -> Result<()> {
    let n_rows = table.nrows();
    let n_cols = table.ncols();
    if contig_lengths.len() != n_rows {
        bail!(
            "Centre log ratio needs one length per row, got {} lengths for {} contigs.",
            contig_lengths.len(),
            n_rows
        );
    }

    let widths = kmer_sizes
        .iter()
        .map(|kmer_size| canonical_count(*kmer_size))
        .collect::<Vec<_>>();
    let wanted = widths.iter().sum::<usize>();
    if wanted != n_cols {
        bail!("k-mer sizes {kmer_sizes:?} want {wanted} columns and the table holds {n_cols}.");
    }

    debug!(
        "Zeros {:.4} of {} composition cells, median replacement {:.3e} at k={}",
        table.iter().filter(|value| **value <= 0.0).count() as f64 / (n_rows * n_cols) as f64,
        n_rows * n_cols,
        median_replacement(contig_lengths, kmer_sizes[0]),
        kmer_sizes[0]
    );

    let transformed = (0..n_rows)
        .into_par_iter()
        .flat_map(|row_index| {
            let row = table.row(row_index);
            let mut centred = Vec::with_capacity(n_cols);
            let mut at = 0;
            for (kmer_size, width) in kmer_sizes.iter().zip(&widths) {
                let block = row.slice(ndarray::s![at..at + width]);
                let block_sum = block.sum();
                let delta = detection_limit(contig_lengths[row_index], *kmer_size);
                let n_zeros = block.iter().filter(|value| **value <= 0.0).count();
                let retained = (1.0 - n_zeros as f64 * delta).max(f64::MIN_POSITIVE);

                let replaced = block
                    .iter()
                    .map(|value| {
                        if *value <= 0.0 || block_sum <= 0.0 {
                            delta
                        } else {
                            value / block_sum * retained
                        }
                    })
                    .collect::<Vec<_>>();

                let log_sum = replaced.iter().map(|value| value.ln()).sum::<f64>();
                let log_geometric_mean = log_sum / *width as f64;
                centred.extend(
                    replaced
                        .into_iter()
                        .map(|value| value.ln() - log_geometric_mean),
                );
                at += width;
            }
            centred
        })
        .collect::<Vec<_>>();

    *table = Array::from_shape_vec((n_rows, n_cols), transformed)?;
    Ok(())
}
