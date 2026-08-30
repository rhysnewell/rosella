//! The two decisions refinement turns on: whether a bin is worth re-clustering at all, and
//! whether the re-clustering that came back is worth taking.

use rosella::clustering::objective::{ClusterObjective, Dbcv};
use rosella::refine::bin_stats::{AGGREGATE, BinStats, EUCLIDEAN, METABAT, RHO, Thresholds};
use rosella::refine::gates::{SplitGate, SplitRejection};
use rosella::refine::splitter::{SplitBars, judge_split, min_validity};

const CONTIG_LENGTH: usize = 100_000;

fn scale() -> rosella::clustering::objective::ScoreThresholds {
    Dbcv::new(&[]).thresholds()
}

fn size_of(cluster: &[usize]) -> usize {
    cluster.len() * CONTIG_LENGTH
}

fn bars(target: f64) -> SplitBars {
    SplitBars {
        target,
        single_cluster: scale().single_cluster,
    }
}

fn cluster(size: usize, offset: usize) -> Vec<usize> {
    (offset..offset + size).collect()
}

#[test]
fn split_rejections() {
    let cases: [(Vec<Vec<usize>>, Vec<usize>, f64, f64, SplitRejection); 5] = [
        (
            vec![cluster(4, 0)],
            vec![],
            1.0,
            0.5,
            SplitRejection::SingleCluster,
        ),
        (
            vec![],
            cluster(4, 0),
            1.0,
            0.5,
            SplitRejection::SingleCluster,
        ),
        (
            vec![cluster(3, 0), cluster(3, 3)],
            vec![],
            0.4,
            0.5,
            SplitRejection::BelowTarget,
        ),
        (
            vec![cluster(3, 0)],
            cluster(1, 3),
            0.5,
            0.0,
            SplitRejection::SingleCluster,
        ),
        (
            vec![cluster(2, 0), cluster(2, 2)],
            cluster(7, 4),
            1.0,
            0.0,
            SplitRejection::AllNoise,
        ),
    ];

    for (clusters, noise, validity, bar, expected) in cases {
        assert_eq!(
            judge_split(
                clusters,
                noise,
                validity,
                bars(bar),
                SplitGate::Strict,
                size_of
            )
            .unwrap_err(),
            expected
        );
    }
}

/// The noise cap is the port's own, not flight's, and it is the one that overrules a
/// re-clustering the density validity already accepted.
#[test]
fn the_validity_gate_takes_a_split_the_noise_cap_rejects() {
    let split = || {
        judge_split(
            vec![cluster(2, 0), cluster(2, 2)],
            cluster(7, 4),
            1.0,
            bars(0.0),
            SplitGate::Validity,
            size_of,
        )
    };

    let (kept, spare) = split().unwrap();
    assert_eq!(kept.len(), 2);
    assert_eq!(spare.len(), 7);
}

/// A piece under the output floor is still a piece. Pouring it in with the noise denies it
/// the recruitment and merge passes that could carry it over the floor.
#[test]
fn a_piece_too_small_to_write_is_still_kept() {
    let (kept, spare) = judge_split(
        vec![cluster(3, 0), cluster(3, 3), cluster(1, 6)],
        cluster(1, 7),
        1.0,
        bars(0.0),
        SplitGate::Strict,
        size_of,
    )
    .unwrap();

    assert_eq!(kept.len(), 3);
    assert_eq!(spare, vec![7]);
}

#[test]
fn a_lone_cluster_survives_on_high_validity() {
    let (kept, spare) = judge_split(
        vec![cluster(3, 0)],
        cluster(1, 3),
        0.95,
        bars(0.0),
        SplitGate::Strict,
        size_of,
    )
    .unwrap();

    assert_eq!(kept.len(), 1);
    assert_eq!(spare, vec![3]);
}

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

const CALM: [f64; 4] = [0.02, 0.02, 1.0, 0.02];
const CALM_THRESHOLDS: [f64; 4] = [0.02, 0.02, 1.0, 0.02];

#[test]
fn a_clean_bin_is_left_alone() {
    let lengths = vec![CONTIG_LENGTH; 20];
    let target = min_validity(
        &stats(CALM, 20),
        &lengths,
        4_000_000,
        false,
        15_000_000,
        &thresholds(CALM_THRESHOLDS),
        scale(),
    );
    assert!(target.is_none());
}

#[test]
fn a_bin_with_too_few_contigs_is_left_alone() {
    let lengths = vec![CONTIG_LENGTH; 9];
    let target = min_validity(
        &stats([0.9, 0.9, 20.0, 0.9], 9),
        &lengths,
        4_000_000,
        false,
        15_000_000,
        &thresholds(CALM_THRESHOLDS),
        scale(),
    );
    assert!(target.is_none());
}

#[test]
fn an_oversized_or_contaminated_bin_splits_on_any_labelling() {
    let lengths = vec![CONTIG_LENGTH; 20];
    let oversized = min_validity(
        &stats(CALM, 20),
        &lengths,
        15_000_000,
        false,
        15_000_000,
        &thresholds(CALM_THRESHOLDS),
        scale(),
    );
    assert_eq!(oversized, Some(0.0));

    let contaminated = min_validity(
        &stats(CALM, 20),
        &lengths,
        4_000_000,
        true,
        15_000_000,
        &thresholds(CALM_THRESHOLDS),
        scale(),
    );
    assert_eq!(contaminated, Some(0.0));
}

/// Both arms of the ladder ask for the same validity. The distance levels that pick between
/// them do not separate a fused bin from a clean one, so a bin's bar cannot rest on them.
#[test]
fn tripped_and_grubby_bins_face_the_same_bar() {
    let lengths = vec![CONTIG_LENGTH; 20];
    let tripped = min_validity(
        &stats([0.5, 0.3, 2.0, 0.5], 20),
        &lengths,
        4_000_000,
        false,
        15_000_000,
        &thresholds(CALM_THRESHOLDS),
        scale(),
    )
    .unwrap();
    let grubby = min_validity(
        &stats([0.1, 0.1, 2.0, 0.2], 20),
        &lengths,
        4_000_000,
        false,
        15_000_000,
        &thresholds(CALM_THRESHOLDS),
        scale(),
    )
    .unwrap();

    assert!(tripped <= 0.5);
    assert_eq!(tripped, grubby);
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

    let lengths = vec![CONTIG_LENGTH; 20];
    let target = min_validity(
        &stats,
        &lengths,
        4_000_000,
        false,
        15_000_000,
        &thresholds(CALM_THRESHOLDS),
        scale(),
    );
    assert!(target.is_some_and(|bar| bar <= 0.5));
}
