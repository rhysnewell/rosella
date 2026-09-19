//! Worth ranks a candidate high enough to claim contigs a cleaner candidate owns. The tiers are
//! a claim order rather than a bar, so the clean candidate is served first and the dirty one is
//! rescored on what is left.

use std::collections::{HashMap, HashSet};

use rosella::clustering::clusterer::Partitioning;
use rosella::clustering::graph_partition::Partition;
use rosella::quality::{Quality, Scorer};
use rosella::recover::ladder::{Judge, combine};
use rosella::recover::peel::peel;

#[path = "../support/bars.rs"]
mod bars;

const CONTIGS: usize = 8;

struct Tabled {
    rows: Vec<(Vec<usize>, f64, f64)>,
}

impl Scorer for Tabled {
    fn score(&self, contigs: &[usize]) -> Quality {
        let mut key = contigs.to_vec();
        key.sort_unstable();
        match self.rows.iter().find(|(members, _, _)| *members == key) {
            Some((_, completeness, contamination)) => Quality {
                completeness: *completeness,
                contamination: *contamination,
                ..Default::default()
            },
            None => Quality::default(),
        }
    }

    fn features(&self, contigs: &[usize]) -> HashSet<u32> {
        contigs.iter().map(|contig| *contig as u32).collect()
    }
}

fn partitioning(clusters: &[&[usize]]) -> Partitioning {
    Partitioning {
        cluster_map: clusters
            .iter()
            .enumerate()
            .map(|(id, contigs)| (id, contigs.iter().copied().collect::<HashSet<_>>()))
            .collect::<HashMap<_, _>>(),
        outliers: HashSet::new(),
        score: None,
        arm: Partition::Leiden,
        seed: 0,
    }
}

fn bins(held: &Partitioning) -> Vec<Vec<usize>> {
    let mut found = held
        .cluster_map
        .values()
        .map(|contigs| {
            let mut contigs = contigs.iter().copied().collect::<Vec<_>>();
            contigs.sort_unstable();
            contigs
        })
        .collect::<Vec<_>>();
    found.sort_unstable();
    found
}

fn arms(clusters: &[&[usize]]) -> Vec<Partitioning> {
    clusters
        .iter()
        .map(|contigs| partitioning(&[contigs]))
        .collect()
}

fn run(
    clusters: &[&[usize]],
    rows: Vec<(Vec<usize>, f64, f64)>,
    lengths: Vec<usize>,
    peeled: bool,
) -> Partitioning {
    let contigs = (0..CONTIGS).collect::<Vec<_>>();
    let scorer = Tabled { rows };
    let judge = Judge {
        quality: &scorer,
        contigs: &contigs,
        bars: bars::bars(80.0),
    };
    match peeled {
        true => peel(&arms(clusters), &judge, &lengths, None),
        false => combine(&arms(clusters), &judge, None),
    }
}

fn even() -> Vec<usize> {
    vec![10_000; CONTIGS]
}

/// Worth ranks the 100/12 candidate over the 70/0 one and lets it take two of the clean
/// candidate's four contigs. The first tier never offers it those contigs at all.
#[test]
fn a_clean_candidate_claims_before_a_dirtier_one_worth_prefers() {
    let clusters: &[&[usize]] = &[&[0, 1, 2, 3], &[2, 3, 4, 5]];
    let rows = vec![
        (vec![0, 1, 2, 3], 100.0, 12.0),
        (vec![2, 3, 4, 5], 70.0, 0.0),
    ];

    assert_eq!(
        bins(&run(clusters, rows.clone(), even(), false)),
        vec![vec![0, 1, 2, 3], vec![4, 5]],
        "worth takes the dirty candidate whole and leaves the clean one a remnant"
    );
    assert_eq!(
        bins(&run(clusters, rows, even(), true)),
        vec![vec![0, 1], vec![2, 3, 4, 5]],
        "the clean candidate is drained first and the dirty one comes back on what is left"
    );
}

/// Two candidates the markers cannot tell apart, so the only thing left to rank them on is how
/// much sequence each one is asking to keep.
#[test]
fn an_equal_ranked_pair_is_broken_by_the_smaller_assembly() {
    let clusters: &[&[usize]] = &[&[0, 1, 2, 3], &[2, 3, 4, 5]];
    let rows = vec![(vec![0, 1, 2, 3], 80.0, 0.0), (vec![2, 3, 4, 5], 80.0, 0.0)];
    let heavy_first = vec![10_000, 10_000, 1_000, 1_000, 500, 500, 500, 500];
    let heavy_second = vec![500, 500, 1_000, 1_000, 10_000, 10_000, 500, 500];

    assert_eq!(
        bins(&run(clusters, rows.clone(), heavy_first, true)),
        vec![vec![0, 1], vec![2, 3, 4, 5]],
        "the second candidate holds less sequence, so it claims the shared contigs"
    );
    assert_eq!(
        bins(&run(clusters, rows, heavy_second, true)),
        vec![vec![0, 1, 2, 3], vec![4, 5]],
        "reversing the sequence reverses the winner, so the break is bases and not order"
    );
}

/// The named tiers stop at 100, and a bin whose markers duplicate over that ceiling is still a
/// bin. Draining only the named tiers would silently unbin it.
#[test]
fn a_candidate_past_the_last_tier_is_still_placed() {
    let held = run(
        &[&[0, 1, 2, 3]],
        vec![(vec![0, 1, 2, 3], 50.0, 150.0)],
        even(),
        true,
    );

    assert_eq!(bins(&held), vec![vec![0, 1, 2, 3]]);
    assert!(
        held.outliers.is_empty(),
        "no tier may drop a candidate the ensemble offered"
    );
}
