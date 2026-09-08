//! Putting bins back in the pot and searching the pool again.

use std::collections::{BTreeMap, HashMap, HashSet};

use ndarray::Array2;
use rosella::clustering::clusterer::HDBSCANResult;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::dissolve::{DissolveScope, DissolveSettings, dissolve};
use rosella::refine::select::Selection;

#[path = "../support/sketches.rs"]
mod sketches;

use sketches::{Fixture, grow, sibling};

const PIECE: usize = 100_000;
const FLOOR: usize = 200_000;
const GENOME: usize = 600_000;

fn settings() -> DissolveSettings {
    DissolveSettings {
        min_bin_size: FLOOR,
        genome_floor: Some(GENOME),
        duplication_bar: 1.0,
        min_contigs: 3,
        scope: DissolveScope::Fused,
        rounds: 1,
        ladder: false,
        n_neighbours: 100,
        completeness: 90.0,
        contamination: 5.0,
        improve: false,
        select: Selection::Rounds,
    }
}

fn result(clusters: Vec<Vec<usize>>, outliers: Vec<usize>) -> HDBSCANResult {
    HDBSCANResult {
        cluster_map: clusters
            .into_iter()
            .enumerate()
            .map(|(id, contigs)| (id, contigs.into_iter().collect::<HashSet<_>>()))
            .collect::<HashMap<_, _>>(),
        outliers: outliers.into_iter().collect(),
        score: 0.0,
    }
}

fn pieces(count: usize) -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    (
        Array2::zeros((count, 2)),
        Array2::zeros((count, 2)),
        vec![PIECE; count],
    )
}

fn bins() -> BTreeMap<usize, Vec<usize>> {
    BTreeMap::from([
        (0usize, vec![0, 1, 2, 3, 4, 5]),
        (1usize, vec![6, 7]),
        (2usize, vec![8]),
    ])
}

/// The pool is only ever committed when a round promotes something, so a run that finds
/// nothing has to leave the bins exactly as the refiner left them.
#[test]
fn a_round_that_promotes_nothing_changes_nothing() {
    let (coverage, tnf, lengths) = pieces(10);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let mut map = bins();
    let mut unbinned = vec![9];

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        settings(),
        |pool, _| {
            Ok(result(
                pool.iter().map(|contig| vec![*contig]).collect(),
                Vec::new(),
            ))
        },
    );

    assert_eq!(ledger.promoted, 0);
    assert_eq!(ledger.refused_small, 4);
    assert_eq!(map, bins());
    assert_eq!(unbinned, vec![9]);
}

/// A bin dissolved into the pool and then not claimed by any accepted cluster has to come back
/// whole. Dropping it strands its long contigs as singletons for nothing.
#[test]
fn a_dissolved_bin_the_pool_does_not_claim_comes_back() {
    let (coverage, tnf, lengths) = pieces(16);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let mut map = BTreeMap::from([
        (0usize, vec![0, 1, 2, 3, 4, 5]),
        (1usize, vec![6, 7]),
        (2usize, vec![8, 9]),
    ]);
    let mut unbinned = vec![10, 11, 12, 13, 14, 15];

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        settings(),
        |_, _| Ok(result(vec![vec![10, 11, 12, 13, 14, 15]], Vec::new())),
    );

    assert_eq!(ledger.promoted, 1);
    assert_eq!(map[&1], vec![6, 7], "an unclaimed bin is restored");
    assert_eq!(map[&2], vec![8, 9]);
    assert!(
        map.values().any(|bin| bin == &vec![10, 11, 12, 13, 14, 15]),
        "the claimed contigs became a bin: {map:?}"
    );
    assert!(unbinned.is_empty());
}

/// One cluster drawing from two bins is the case where a contig can end up in the new bin and
/// in its old one at once, and where the ledger has to close: everything the pool held is
/// adopted, returned or left.
#[test]
fn a_cluster_drawing_from_two_bins_leaves_neither_holding_it() {
    let (coverage, tnf, lengths) = pieces(16);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let mut map = BTreeMap::from([
        (0usize, vec![0, 1, 2, 3, 4, 5]),
        (1usize, vec![6, 7, 8, 9]),
        (2usize, vec![10, 11, 12, 13]),
    ]);
    let mut unbinned = vec![14, 15];

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        settings(),
        |_, _| Ok(result(vec![vec![6, 7, 10, 11, 14, 15]], Vec::new())),
    );

    assert_eq!(
        map[&0],
        vec![0, 1, 2, 3, 4, 5],
        "a whole bin never dissolves"
    );
    assert_eq!(map[&1], vec![8, 9]);
    assert_eq!(map[&2], vec![12, 13]);
    assert!(map.values().any(|bin| bin == &vec![6, 7, 10, 11, 14, 15]));
    assert!(unbinned.is_empty());

    assert_eq!(ledger.pool_contigs, 10);
    assert_eq!(
        ledger.adopted_contigs + ledger.returned_contigs + ledger.left_contigs,
        ledger.pool_contigs,
        "every contig the pool held has one fate: {ledger}"
    );
    assert_eq!(
        ledger.adopted_bp + ledger.returned_bp + ledger.left_bp,
        ledger.pool_bp
    );
}

/// The pool's partition is inserted as bins with nothing else checking it, so a contig in two
/// clusters would reach the writer, where one label silently wins.
#[test]
fn a_pool_that_places_a_contig_twice_is_refused() {
    let (coverage, tnf, lengths) = pieces(16);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let start = BTreeMap::from([(0usize, vec![0, 1, 2, 3, 4, 5]), (1usize, vec![6, 7, 8, 9])]);
    let mut map = start.clone();
    let mut unbinned = vec![10, 11, 12, 13, 14, 15];

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        settings(),
        |_, _| {
            Ok(result(
                vec![vec![6, 7, 8, 9, 10, 11], vec![10, 11, 12, 13, 14, 15]],
                Vec::new(),
            ))
        },
    );

    assert_eq!(ledger.proposed, 2);
    assert_eq!(ledger.promoted, 0);
    assert_eq!(map, start);
    assert_eq!(unbinned, vec![10, 11, 12, 13, 14, 15]);
}

fn fused() -> Fixture {
    let native = grow(3, 300_000);
    Fixture::new(
        "dissolve_fused",
        vec![
            native.clone(),
            sibling(&native[0..30_000]),
            sibling(&native[100_000..130_000]),
            sibling(&native[200_000..230_000]),
            grow(41, 50_000),
        ],
    )
}

/// The eject arm has already had its go, so a bin still holding its own sequence twice goes in
/// the pot on those grounds alone, and the same sequence cannot buy its way back out.
#[test]
fn a_bin_holding_its_own_sequence_twice_dissolves_and_is_not_taken_back() {
    let fixture = fused();
    let features = fixture.features();
    let start = BTreeMap::from([(0usize, vec![0, 1, 2, 3])]);
    let mut map = start.clone();
    let mut unbinned = vec![4];
    let settings = DissolveSettings {
        genome_floor: Some(300_000),
        duplication_bar: 0.05,
        ..settings()
    };

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        settings,
        |_, _| Ok(result(vec![vec![0, 1, 2, 3], vec![4]], Vec::new())),
    );

    assert_eq!(ledger.dissolved_duplicated, 1, "{ledger}");
    assert_eq!(ledger.dissolved_small, 0);
    assert_eq!(ledger.refused_duplicated, 1);
    assert_eq!(ledger.refused_small, 1);
    assert_eq!(ledger.promoted, 0);
    assert_eq!(map, start);
}

/// Each round searches only what the round before it left, at half the neighbours, which is the
/// whole point of the loop: the graph is rebuilt over material that changed.
#[test]
fn a_round_searches_what_the_round_before_it_left() {
    let (coverage, tnf, lengths) = pieces(16);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let mut map = BTreeMap::from([(0usize, vec![0, 1]), (1usize, vec![2, 3])]);
    let mut unbinned = (4..16).collect::<Vec<_>>();
    let seen = std::cell::RefCell::new(Vec::new());

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        DissolveSettings {
            rounds: 3,
            ..settings()
        },
        |pool, round| {
            let mut pool = pool.iter().copied().collect::<Vec<_>>();
            pool.sort_unstable();
            seen.borrow_mut().push((pool.clone(), round.n_neighbours));
            let take = pool.iter().copied().take(6).collect::<Vec<_>>();
            Ok(result(vec![take], Vec::new()))
        },
    );

    let seen = seen.into_inner();
    assert_eq!(ledger.rounds, 3, "{ledger}");
    assert_eq!(
        ledger.promoted, 2,
        "the last round is left 4 contigs, under the floor"
    );
    assert_eq!(
        seen.iter().map(|(_, k)| *k).collect::<Vec<_>>(),
        vec![100, 50, 25],
        "the neighbour count halves as the pool shrinks"
    );
    assert_eq!(seen[0].0.len(), 16);
    assert_eq!(seen[1].0.len(), 10);
    assert_eq!(seen[2].0.len(), 4);
    assert!(
        seen[1]
            .0
            .iter()
            .all(|contig| !seen[0].0[..6].contains(contig)),
        "an accepted cluster leaves the pool"
    );
}

/// `all` puts a clean bin over the genome floor in the pot too, and the return rule is the only
/// thing that makes that safe: a bin no cluster claims has to come back whole.
#[test]
fn scope_all_dissolves_a_clean_bin_and_gives_it_back_unclaimed() {
    let (coverage, tnf, lengths) = pieces(16);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let clean = (0..8).collect::<Vec<_>>();
    let mut map = BTreeMap::from([(0usize, clean.clone()), (1usize, vec![8, 9])]);
    let mut unbinned = (10..16).collect::<Vec<_>>();

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        DissolveSettings {
            scope: DissolveScope::All,
            ..settings()
        },
        |_, _| Ok(result(vec![vec![10, 11, 12, 13, 14, 15]], Vec::new())),
    );

    assert_eq!(ledger.dissolved_clean, 1, "{ledger}");
    assert_eq!(ledger.pool_contigs, 16);
    assert_eq!(map[&0], clean, "the clean bin came back whole");
    assert_eq!(map[&1], vec![8, 9]);
}

/// A round that accepts nothing at the genome floor is where the loop stops, so the ladder is
/// what decides whether it stops there or drops the floor and takes the same clusters.
#[test]
fn the_ladder_takes_clusters_the_fixed_floor_refuses() {
    let (coverage, tnf, lengths) = pieces(16);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let start = BTreeMap::from([(0usize, vec![0, 1]), (1usize, vec![2, 3])]);
    let clusters = vec![(4..9).collect::<Vec<_>>(), (9..14).collect::<Vec<_>>()];

    let run = |ladder: bool| {
        let mut map = start.clone();
        let mut unbinned = (4..16).collect::<Vec<_>>();
        let ledger = dissolve(
            &features,
            None,
            &mut map,
            &mut unbinned,
            DissolveSettings {
                genome_floor: Some(GENOME),
                ladder,
                ..settings()
            },
            |_, _| Ok(result(clusters.clone(), Vec::new())),
        );
        (ledger, map)
    };

    let (fixed, _) = run(false);
    assert_eq!(
        fixed.promoted, 0,
        "500 kb clusters under a 600 kb genome floor"
    );
    assert_eq!(fixed.refused_small, 2);

    let (relaxed, map) = run(true);
    assert_eq!(relaxed.promoted, 2, "{relaxed}");
    assert!(relaxed.rung > 0);
    assert!(map.values().any(|bin| bin == &clusters[0]));
}

/// Ranked selection lets every round propose over the same pool, so the rule that decides which
/// proposal survives is the bar and the claim, not which round ran first.
#[test]
fn ranked_selection_gives_a_contig_to_one_proposal_only() {
    let (coverage, tnf, lengths) = pieces(16);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let mut map = BTreeMap::from([(0usize, vec![0, 1, 2, 3, 4, 5, 6, 7])]);
    let mut unbinned = (8..16).collect::<Vec<_>>();
    let settings = DissolveSettings {
        scope: DissolveScope::All,
        rounds: 2,
        select: Selection::Ranked,
        ..settings()
    };

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        settings,
        |_, round| {
            Ok(match round.n_neighbours {
                100 => result(vec![(0..12).collect()], Vec::new()),
                _ => result(vec![(0..6).collect(), (6..12).collect()], Vec::new()),
            })
        },
    );

    assert_eq!(ledger.rounds, 2);
    assert_eq!(ledger.proposed, 3);
    assert_eq!(
        ledger.promoted, 1,
        "the halves overlap what the blob claimed"
    );
    let placed = map.values().flatten().copied().collect::<Vec<_>>();
    let mut seen = placed.clone();
    seen.sort_unstable();
    seen.dedup();
    assert_eq!(seen.len(), placed.len(), "no contig lands in two bins");
}

/// A proposal the winner did not touch still stands, so ranking costs nothing where the rounds
/// do not disagree.
#[test]
fn ranked_selection_keeps_proposals_that_do_not_overlap() {
    let (coverage, tnf, lengths) = pieces(20);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let mut map = BTreeMap::from([(0usize, vec![0, 1, 2, 3, 4, 5, 6, 7])]);
    let mut unbinned = (8..20).collect::<Vec<_>>();
    let settings = DissolveSettings {
        scope: DissolveScope::All,
        rounds: 2,
        select: Selection::Ranked,
        ..settings()
    };

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        settings,
        |_, round| {
            Ok(match round.n_neighbours {
                100 => result(vec![(0..8).collect()], Vec::new()),
                _ => result(vec![(8..16).collect()], Vec::new()),
            })
        },
    );

    assert_eq!(ledger.promoted, 2);
    assert_eq!(ledger.adopted_contigs, 16);
}
