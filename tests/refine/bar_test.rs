//! Whether a bin is worth re-clustering at all.

use rosella::refine::bar::should_split;
use rosella::refine::bin_stats::{ABUNDANCE, AGGREGATE, BinStats, EUCLIDEAN, RHO, Thresholds};
use rosella::refine::gates::Trigger;

const CONTIG_LENGTH: usize = 100_000;
const MAX_BIN_SIZE: usize = 15_000_000;

const CALM: [f64; 4] = [0.02, 0.02, 1.0, 0.02];
const CALM_THRESHOLDS: [f64; 4] = [0.05, 0.05, 1.1, 0.05];

fn stats(mean: [f64; 4], n: usize) -> BinStats {
    BinStats {
        mean,
        std: [0.1; 4],
        per_contig: vec![mean; n],
    }
}

fn thresholds(quantile: [f64; 4]) -> Thresholds {
    Thresholds { quantile }
}

fn bar(stats: &BinStats, contigs: usize, bin_size: usize) -> Option<Trigger> {
    over_budget(stats, contigs, bin_size, false)
}

fn over_budget(
    stats: &BinStats,
    contigs: usize,
    bin_size: usize,
    over_budget: bool,
) -> Option<Trigger> {
    should_split(
        stats,
        &vec![CONTIG_LENGTH; contigs],
        bin_size,
        over_budget,
        MAX_BIN_SIZE,
        &thresholds(CALM_THRESHOLDS),
    )
}

#[test]
fn a_clean_bin_is_left_alone() {
    assert!(bar(&stats(CALM, 20), 20, 4_000_000).is_none());
}

#[test]
fn a_bin_with_too_few_contigs_is_left_alone() {
    assert!(bar(&stats([0.9, 0.9, 20.0, 0.9], 9), 9, 4_000_000).is_none());
}

#[test]
fn an_oversized_or_contaminated_bin_splits_on_any_labelling() {
    let oversized = bar(&stats(CALM, 20), 20, MAX_BIN_SIZE);
    assert_eq!(oversized, Some(Trigger::Forced));

    let contaminated = over_budget(&stats(CALM, 20), 20, 4_000_000, true);
    assert_eq!(contaminated, Some(Trigger::Forced));
}

/// Bin averages can sit under every level while individual contigs sit well over them.
/// Enough of those by length is a trigger on its own.
#[test]
fn misplaced_contigs_trigger_on_their_own() {
    let mut stats = stats(CALM, 20);
    for row in stats.per_contig.iter_mut().take(11) {
        row[ABUNDANCE] = 0.9;
        row[RHO] = 0.9;
        row[EUCLIDEAN] = 30.0;
        row[AGGREGATE] = 0.9;
    }

    assert_eq!(
        bar(&stats, 20, 4_000_000).unwrap(),
        Trigger::Tripped {
            columns: [false; 4],
            misplaced: true
        }
    );
}

#[test]
fn only_coverage_keeps_a_floor() {
    let mut learned = CALM;
    learned[EUCLIDEAN] = 1.2;
    assert_eq!(
        bar(&stats(learned, 20), 20, 4_000_000),
        Some(Trigger::Tripped {
            columns: [false, false, true, false],
            misplaced: true
        })
    );

    let mut depth = CALM;
    depth[ABUNDANCE] = 0.2;
    assert!(bar(&stats(depth, 20), 20, 4_000_000).is_none());
}
