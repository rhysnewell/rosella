use crate::refine::bin_stats::{AGGREGATE, BinStats, EUCLIDEAN, METABAT, RHO, Thresholds};
use crate::refine::gates::Trigger;

/// Floors on each level, so a run where every bin looks alike does not start splitting on
/// noise. flight's `validate_bins`.
const FLOORS: [f64; 4] = [0.30, 0.15, 6.0, 0.35];
const MULTIPLIERS: [f64; 4] = [1.25, 1.5, 1.25, 1.5];

/// Contigs flagged as out of place have to add up to this before they alone trigger a
/// split.
const MISPLACED_LENGTH: usize = 1_000_000;

/// Below this there is not enough of a bin to re-cluster, so the work is wasted.
pub const MIN_SPLIT_CONTIGS: usize = 10;

/// flight reads its rungs off absolute genome sizes against a 20 Mbp ceiling. Held as
/// fractions instead so the ladder follows `--max-bin-size` rather than assuming the
/// assembly holds bacteria of a particular size.
const RUNGS: [(f64, Option<usize>, f64); 3] = [
    (0.8, None, 2.5),
    (0.7, Some(2 * MISPLACED_LENGTH), 2.0),
    (0.6, Some(MISPLACED_LENGTH), 1.5),
];

/// Only bins that trip a level are re-clustered. Merely dirty bins were 62% of candidates and
/// 18% of accepted splits, and dropping them raised t1 on CAMI I high.
pub fn min_validity(
    stats: &BinStats,
    lengths: &[usize],
    bin_size: usize,
    over_budget: bool,
    max_bin_size: usize,
    thresholds: &Thresholds,
) -> Option<(f64, Trigger)> {
    if lengths.len() < MIN_SPLIT_CONTIGS {
        return None;
    }
    if bin_size >= max_bin_size || over_budget {
        return Some((0.0, Trigger::Forced));
    }

    let levels = levels(thresholds);
    let mut columns = [false; 4];
    for (column, fired) in columns.iter_mut().enumerate() {
        *fired = stats.mean[column] >= levels[column];
    }
    let misplaced = misplaced_length(stats, lengths, &levels);

    if columns.iter().any(|fired| *fired) || misplaced >= MISPLACED_LENGTH {
        let ceiling = levels[METABAT].max(levels[AGGREGATE]);
        let factor = match rung(bin_size, max_bin_size, misplaced) {
            Some(multiplier) => ceiling * multiplier,
            None => spread(stats) * 1.25,
        };
        return Some((
            (1.0 - factor.min(1.0)).clamp(0.0, 1.0),
            Trigger::Tripped {
                columns,
                misplaced: misplaced >= MISPLACED_LENGTH,
            },
        ));
    }

    None
}

/// Rungs read the cross-bin levels, the fallthrough reads the bin's own spread: landing on a
/// rung says the bin is big or ragged, neither of which says how tight it is.
fn rung(bin_size: usize, max_bin_size: usize, misplaced: usize) -> Option<f64> {
    let fraction = bin_size as f64 / max_bin_size as f64;
    RUNGS
        .iter()
        .find(|(size, length, _)| {
            fraction >= *size || length.is_some_and(|floor| misplaced >= floor)
        })
        .map(|(_, _, multiplier)| *multiplier)
}

fn spread(stats: &BinStats) -> f64 {
    (stats.mean[AGGREGATE] + stats.std[AGGREGATE]).max(stats.mean[METABAT] + stats.std[METABAT])
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

/// flight caps each level at `mean + std * 1.5`. That cap cannot change the mean test, but it
/// does bind in `misplaced_length`, where flight compares single contigs against the level.
fn levels(thresholds: &Thresholds) -> [f64; 4] {
    let mut levels = [0.0f64; 4];
    for (column, level) in levels.iter_mut().enumerate() {
        *level = FLOORS[column].max(MULTIPLIERS[column] * thresholds.mean[column]);
    }
    levels
}
