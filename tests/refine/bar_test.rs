//! Whether a bin is worth re-clustering at all, and what a re-clustering of it has to reach.

use rosella::refine::bar::min_validity;
use rosella::refine::bin_stats::{AGGREGATE, BinStats, EUCLIDEAN, METABAT, RHO, Thresholds};
use rosella::refine::gates::Trigger;

const CONTIG_LENGTH: usize = 100_000;
const MAX_BIN_SIZE: usize = 15_000_000;

const CALM: [f64; 4] = [0.02, 0.02, 1.0, 0.02];
const CALM_THRESHOLDS: [f64; 4] = [0.02, 0.02, 1.0, 0.02];

fn stats(mean: [f64; 4], n: usize) -> BinStats {
    BinStats {
        mean,
        std: [0.1; 4],
        per_contig: vec![mean; n],
    }
}

fn thresholds(mean: [f64; 4]) -> Thresholds {
    Thresholds { mean }
}

fn bar(stats: &BinStats, contigs: usize, bin_size: usize) -> Option<(f64, Trigger)> {
    over_budget(stats, contigs, bin_size, false)
}

fn over_budget(
    stats: &BinStats,
    contigs: usize,
    bin_size: usize,
    over_budget: bool,
) -> Option<(f64, Trigger)> {
    min_validity(
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
    assert_eq!(oversized, Some((0.0, Trigger::Forced)));

    let contaminated = over_budget(&stats(CALM, 20), 20, 4_000_000, true);
    assert_eq!(contaminated, Some((0.0, Trigger::Forced)));
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

    let (target, trigger) = bar(&stats, 20, 4_000_000).unwrap();
    assert!(target <= 0.5);
    assert_eq!(
        trigger,
        Trigger::Tripped {
            columns: [false; 4],
            misplaced: true
        }
    );
}

/// Tripped on the aggregate mean alone, which `misplaced_length` does not read, so the rung
/// is picked by size and nothing else.
fn tripped_on_aggregate(std: f64) -> BinStats {
    BinStats {
        mean: [0.1, 0.1, 2.0, 0.4],
        std: [std; 4],
        per_contig: vec![[0.1, 0.1, 2.0, 0.4]; 20],
    }
}

#[test]
fn the_bar_falls_as_the_bin_grows_across_the_rungs() {
    let stats = tripped_on_aggregate(0.1);
    let rungs = [9_000_000, 10_500_000, 12_000_000, MAX_BIN_SIZE];

    let bars = rungs
        .iter()
        .map(|size| bar(&stats, 20, *size).unwrap().0)
        .collect::<Vec<_>>();

    assert!(
        bars.windows(2).all(|pair| pair[0] > pair[1]),
        "rungs did not fall: {bars:?}"
    );
}

/// A rung is a constant off the cross-bin levels while the fallthrough reads the bin, so
/// widening the bin moves one bar and not the other.
#[test]
fn only_the_fallthrough_follows_the_bins_own_spread() {
    let tight = tripped_on_aggregate(0.05);
    let wide = tripped_on_aggregate(0.3);

    let under_every_rung = 8_000_000;
    let on_a_rung = 9_000_000;

    assert_ne!(
        bar(&tight, 20, under_every_rung).unwrap().0,
        bar(&wide, 20, under_every_rung).unwrap().0
    );
    assert_eq!(
        bar(&tight, 20, on_a_rung).unwrap().0,
        bar(&wide, 20, on_a_rung).unwrap().0
    );
}

/// The two formulas are not comparable, so the ladder inverts where they meet: a ragged bin
/// one rung short faces a lower bar than the same bin on the rung. flight's, and it needs a
/// spread over 0.42 to show up at all.
#[test]
fn the_ladder_inverts_where_the_fallthrough_meets_the_first_rung() {
    let ragged = tripped_on_aggregate(0.1);

    let under = bar(&ragged, 20, 8_000_000).unwrap().0;
    let on_rung = bar(&ragged, 20, 9_000_000).unwrap().0;
    assert!(under < on_rung, "{under} should undercut {on_rung}");

    let tight = tripped_on_aggregate(0.0);
    assert!(bar(&tight, 20, 8_000_000).unwrap().0 > bar(&tight, 20, 9_000_000).unwrap().0);
}
