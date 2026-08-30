//! Whether the re-clustering that came back is worth taking. Whether the bin was worth
//! re-clustering in the first place is `bar_test`.

use rosella::clustering::objective::{ClusterObjective, Dbcv};
use rosella::refine::gates::{SplitGate, SplitRejection};
use rosella::refine::splitter::{SplitBars, judge_split};

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
