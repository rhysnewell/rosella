use crate::refine::bin_stats::BinStats;
use crate::refine::gates::{SplitRejection, Trigger};

pub(crate) enum Proposal {
    NoStats,
    TooFewContigs,
    NoTrigger,
    NoClustering(Trigger),
    Rejected(Trigger, SplitRejection),
    Accepted(Trigger, SplitOutcome),
}

pub(crate) struct SplitOutcome {
    pub kept: Vec<Vec<usize>>,
    pub unbinned: Vec<usize>,
}

/// Size is not a bar here. A piece too small to write out can still recruit or merge its way
/// over the floor, so `bin_writer` applies `min_bin_size` once, at the end.
pub fn judge_split(
    clusters: Vec<Vec<usize>>,
    noise: Vec<usize>,
    size_of: impl Fn(&[usize]) -> usize,
) -> Result<(Vec<Vec<usize>>, Vec<usize>), SplitRejection> {
    let distinct = clusters.len() + usize::from(!noise.is_empty());
    if distinct <= 1 {
        return Err(SplitRejection::SingleCluster);
    }

    let bin_size = clusters
        .iter()
        .chain(std::iter::once(&noise))
        .map(|contigs| size_of(contigs))
        .sum::<usize>() as f64;
    if size_of(&noise) as f64 > crate::tuning::MAX_NOISE_FRACTION * bin_size {
        return Err(SplitRejection::AllNoise);
    }

    Ok((clusters, noise))
}

pub fn leaves_two_standing(
    pieces: &[Vec<usize>],
    floor: usize,
    size_of: impl Fn(&[usize]) -> usize,
) -> bool {
    pieces
        .iter()
        .filter(|piece| size_of(piece) >= floor)
        .count()
        >= 2
}

pub(crate) fn tighter(pieces: f64, whole: &BinStats, column: usize) -> bool {
    pieces <= whole.mean[column] * crate::tuning::REQUIRED_IMPROVEMENT
}

pub(crate) fn contigs(indices: &[usize], positions: impl Iterator<Item = usize>) -> Vec<usize> {
    let mut contigs = positions
        .map(|position| indices[position])
        .collect::<Vec<_>>();
    contigs.sort_unstable();
    contigs
}
