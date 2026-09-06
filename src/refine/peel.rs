use crate::refine::bin_stats::{AGGREGATE, BinStats};

pub struct Peel {
    pub lone: Vec<usize>,
    pub rest: Vec<usize>,
}

/// A contig long enough to report alone that sits further from its bin than the bin's own
/// spread is the closed genome a partition absorbs, and re-clustering the bin at the same k
/// only surrounds it with the same neighbours again.
pub fn candidate(
    indices: &[usize],
    stats: &BinStats,
    lengths: &[usize],
    min_bin_size: usize,
) -> Option<Peel> {
    let bar = stats.mean[AGGREGATE] + stats.std[AGGREGATE];
    let mut lone = Vec::new();
    let mut rest = Vec::new();
    for ((contig, length), row) in indices.iter().zip(lengths).zip(&stats.per_contig) {
        if *length >= min_bin_size && row[AGGREGATE] > bar {
            lone.push(*contig);
        } else {
            rest.push(*contig);
        }
    }
    if lone.is_empty() || rest.len() < 2 {
        return None;
    }
    Some(Peel { lone, rest })
}
