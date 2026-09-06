//! Whether a bin is worth re-clustering at all.

use rosella::refine::bar::{describe_levels, should_split};
use rosella::refine::bin_stats::{
    AGGREGATE, BinStats, EUCLIDEAN, LevelSource, METABAT, RHO, Thresholds,
};
use rosella::refine::gates::Trigger;

const CONTIG_LENGTH: usize = 100_000;
const MAX_BIN_SIZE: usize = 15_000_000;

const CALM: [f64; 4] = [0.02, 0.02, 1.0, 0.02];
const REALISTIC: [f64; 4] = [0.072, 0.045, 2.1, 0.068];
const CALM_THRESHOLDS: [f64; 4] = [0.02, 0.02, 1.0, 0.02];

fn stats(mean: [f64; 4], n: usize) -> BinStats {
    BinStats {
        mean,
        std: [0.1; 4],
        per_contig: vec![mean; n],
    }
}

fn thresholds(mean: [f64; 4]) -> Thresholds {
    Thresholds {
        mean,
        source: LevelSource::Flight,
    }
}

fn derived(mean: [f64; 4]) -> Thresholds {
    Thresholds {
        mean,
        source: LevelSource::Derived,
    }
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
        row[METABAT] = 0.9;
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

/// The recorded bin spreads across the three CAMI sets are 0.072, 0.094 and 0.114, so under
/// flight's constants every level is a floor and none of them describes the run.
#[test]
fn a_realistic_run_never_reaches_a_level_off_its_own_means() {
    let described = describe_levels(&thresholds(REALISTIC));
    assert_eq!(described.matches("(floor,").count(), 4, "{described}");
    assert!(
        described.contains("aggregate 0.3500 (floor, mean 0.0680)"),
        "{described}"
    );
}

/// 0.35 is a geometric-scale number and the aggregate column is the only one that passes
/// through the combination, so it is the only level the arithmetic default has to move.
#[test]
fn derived_levels_put_the_aggregate_column_where_the_arithmetic_distance_lives() {
    let derived = describe_levels(&derived(REALISTIC));
    assert!(derived.contains("aggregate 0.1225 (floor,"), "{derived}");
    assert!(derived.contains("metabat 0.3000 (floor,"), "{derived}");

    let clears = describe_levels(&self::derived([0.31, 0.16, 6.1, 0.14]));
    assert_eq!(clears.matches("(run,").count(), 4, "{clears}");
}
