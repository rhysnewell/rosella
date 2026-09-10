use crate::refine::bin_stats::{BinStats, EUCLIDEAN, LevelSource, METABAT, RHO, Thresholds};
use crate::refine::gates::Trigger;

/// Floors on each level, so a run where every bin looks alike does not start splitting on
/// noise. flight's `validate_bins`.
const FLOORS: [f64; 4] = [0.30, 0.15, 6.0, 0.35];
const MULTIPLIERS: [f64; 4] = [1.25, 1.5, 1.25, 1.5];

/// Only the aggregate column passes through `Combination::combine`, whose geometric branch
/// carries a square root the arithmetic default does not, so flight's 0.35 is a 0.1225 here.
const GUARDS: [f64; 4] = [0.30, 0.15, 6.0, 0.1225];

/// Contigs flagged as out of place have to add up to this before they alone trigger a
/// split.
const MISPLACED_LENGTH: usize = 1_000_000;

/// Below this there is not enough of a bin to re-cluster, so the work is wasted.
pub const MIN_SPLIT_CONTIGS: usize = 10;

/// Only bins that trip a level are re-clustered. Merely dirty bins were 62% of candidates and
/// 18% of accepted splits, and dropping them raised t1 on CAMI I high.
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

    if columns.iter().any(|fired| *fired) || misplaced >= MISPLACED_LENGTH {
        return Some(Trigger::Tripped {
            columns,
            misplaced: misplaced >= MISPLACED_LENGTH,
        });
    }

    None
}

fn misplaced_length(stats: &BinStats, lengths: &[usize], levels: &[f64; 4]) -> usize {
    lengths
        .iter()
        .zip(stats.per_contig.iter())
        .filter(|(_, averages)| {
            averages[METABAT] >= levels[METABAT]
                || averages[RHO] >= levels[RHO]
                || averages[EUCLIDEAN] >= levels[EUCLIDEAN]
        })
        .map(|(length, _)| *length)
        .sum()
}

const COLUMN_NAMES: [&str; 4] = ["metabat", "rho", "euclidean", "aggregate"];

/// The floors came from flight's geometric distances, so on an arithmetic run they can pin a
/// column open or shut, and no log said whether the floor or the run's own mean bound.
pub fn describe_levels(thresholds: &Thresholds) -> String {
    let levels = levels(thresholds);
    let floors = floors(thresholds.source);
    COLUMN_NAMES
        .iter()
        .enumerate()
        .map(|(column, name)| {
            let source = if levels[column] > floors[column] {
                "run"
            } else {
                "floor"
            };
            format!(
                "{name} {:.4} ({source}, mean {:.4})",
                levels[column], thresholds.mean[column]
            )
        })
        .collect::<Vec<_>>()
        .join(", ")
}

/// flight's multiplier lifts a mean to "worse than typical". A quantile already is that, so
/// applying one on top of `Derived` would count the same allowance twice.
pub fn levels(thresholds: &Thresholds) -> [f64; 4] {
    let floors = floors(thresholds.source);
    let mut levels = [0.0f64; 4];
    for (column, level) in levels.iter_mut().enumerate() {
        *level = match thresholds.source {
            LevelSource::Flight => {
                floors[column].max(MULTIPLIERS[column] * thresholds.mean[column])
            }
            LevelSource::Derived => floors[column].max(thresholds.mean[column]),
        };
    }
    levels
}

fn floors(source: LevelSource) -> [f64; 4] {
    match source {
        LevelSource::Flight => FLOORS,
        LevelSource::Derived => GUARDS,
    }
}
