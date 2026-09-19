use anyhow::{Result, bail};
use log::debug;
use ndarray::{Array, Array2};
use rayon::prelude::*;

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

    debug!(
        "Zeros {:.4} of {} composition cells, median replacement {:.3e} at k={}",
        table.iter().filter(|value| **value <= 0.0).count() as f64 / (n_rows * n_cols) as f64,
        n_rows * n_cols,
        median_replacement(contig_lengths, kmer_size),
        kmer_size
    );

    let transformed = (0..n_rows)
        .into_par_iter()
        .flat_map(|row_index| {
            let row = table.row(row_index);
            let row_sum = row.sum();
            let delta = detection_limit(contig_lengths[row_index], kmer_size);
            let n_zeros = row.iter().filter(|value| **value <= 0.0).count();
            let retained = (1.0 - n_zeros as f64 * delta).max(f64::MIN_POSITIVE);

            let replaced = row
                .iter()
                .map(|value| {
                    if *value <= 0.0 || row_sum <= 0.0 {
                        delta
                    } else {
                        value / row_sum * retained
                    }
                })
                .collect::<Vec<_>>();

            let log_geometric_mean =
                replaced.iter().map(|value| value.ln()).sum::<f64>() / n_cols as f64;
            replaced
                .into_iter()
                .map(|value| value.ln() - log_geometric_mean)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    *table = Array::from_shape_vec((n_rows, n_cols), transformed)?;
    Ok(())
}
