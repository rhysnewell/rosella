use crate::refine::bin_stats::{ABUNDANCE, BinStats, EUCLIDEAN, RHO, Thresholds};
use crate::refine::gates::Trigger;

// Only coverage kept its scale when the learned rows replaced the k-mers, so only it keeps a floor.
const ABUNDANCE_FLOOR: f64 = 0.30;

// Below this there is not enough of a bin to re-cluster, so the work is wasted. Swept over 14
// AMBER sets: 3 and 5 are byte identical to 10, and 20 loses t1, so nothing in the range decides.
pub const MIN_SPLIT_CONTIGS: usize = 10;

// Only bins that trip a level are re-clustered. Merely dirty bins were 62% of candidates and
// 18% of accepted splits, and dropping them raised t1 on CAMI I high.
pub fn should_split(
    stats: &BinStats,
    lengths: &[usize],
    bin_size: usize,
    over_budget: bool,
    max_bin_size: usize,
    thresholds: &Thresholds,
) -> Option<Trigger> {
    if lengths.len() < MIN_SPLIT_CONTIGS {
        return None;
    }
    if bin_size >= max_bin_size || over_budget {
        return Some(Trigger::Forced);
    }

    let levels = levels(thresholds);
    let mut columns = [false; 4];
    for (column, fired) in columns.iter_mut().enumerate() {
        *fired = stats.mean[column] >= levels[column];
    }
    let misplaced = misplaced_length(stats, lengths, &levels);

    if columns.iter().any(|fired| *fired) || misplaced >= crate::tuning::MISPLACED_LENGTH {
        return Some(Trigger::Tripped {
            columns,
            misplaced: misplaced >= crate::tuning::MISPLACED_LENGTH,
        });
    }

    None
}

fn misplaced_length(stats: &BinStats, lengths: &[usize], levels: &[f64; 4]) -> usize {
    lengths
        .iter()
        .zip(stats.per_contig.iter())
        .filter(|(_, averages)| {
            averages[ABUNDANCE] >= levels[ABUNDANCE]
                || averages[RHO] >= levels[RHO]
                || averages[EUCLIDEAN] >= levels[EUCLIDEAN]
        })
        .map(|(length, _)| *length)
        .sum()
}

const COLUMN_NAMES: [&str; 4] = ["abundance", "rho", "euclidean", "aggregate"];

pub fn describe_levels(thresholds: &Thresholds) -> String {
    let levels = levels(thresholds);
    COLUMN_NAMES
        .iter()
        .zip(levels)
        .map(|(name, level)| format!("{name} {level:.4}"))
        .collect::<Vec<_>>()
        .join(", ")
}

fn levels(thresholds: &Thresholds) -> [f64; 4] {
    let mut levels = thresholds.quantile;
    levels[ABUNDANCE] = levels[ABUNDANCE].max(ABUNDANCE_FLOOR);
    levels
}
