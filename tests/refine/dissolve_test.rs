//! Putting bins back in the pot and searching the pool again.

use std::collections::{BTreeMap, HashMap, HashSet};

use ndarray::Array2;
use rosella::clustering::clusterer::Partitioning;
use rosella::embedding::features::ContigFeatures;
use rosella::embedding::knn::KnnGraph;
use rosella::refine::dissolve::{DissolveSettings, POOL_VIEWS, dissolve};
use rosella::refine::rung::Bars;

#[path = "../support/sketches.rs"]
mod sketches;

use sketches::{Fixture, grow, sibling};

const PIECE: usize = 100_000;
const FLOOR: usize = 200_000;
const GENOME: usize = 600_000;

fn settings() -> DissolveSettings {
    DissolveSettings {
        bars: Bars {
            min_bin_size: FLOOR,
            duplication_bar: 1.0,
            completeness: 90.0,
            contamination: 5.0,
        },
        genome_floor: Some(GENOME),
        min_contigs: 3,
        rounds: 1,
        passes: 1,
        n_neighbours: NEIGHBOURS,
        reuse: true,
    }
}

fn partitioning(clusters: Vec<Vec<usize>>, outliers: Vec<usize>) -> Partitioning {
    Partitioning {
        cluster_map: clusters
            .into_iter()
            .enumerate()
            .map(|(id, contigs)| (id, contigs.into_iter().collect::<HashSet<_>>()))
            .collect::<HashMap<_, _>>(),
        outliers: outliers.into_iter().collect(),
        score: 0.0,
    }
}

fn result(clusters: Vec<Vec<usize>>, outliers: Vec<usize>) -> Vec<Partitioning> {
    vec![partitioning(clusters, outliers)]
}

/// The rounds truncate what the search hands back, so a pinned search still needs a real graph.
const NEIGHBOURS: usize = 100;

fn empty_knn(pool: &HashSet<usize>) -> Result<(KnnGraph, Vec<usize>), anyhow::Error> {
    let mut order = pool.iter().copied().collect::<Vec<_>>();
    order.sort_unstable();
    Ok((
        KnnGraph {
            indices: Array2::zeros((order.len(), NEIGHBOURS)),
            dists: Array2::zeros((order.len(), NEIGHBOURS)),
        },
        order,
    ))
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
        &[],
        |pool, _, _| empty_knn(pool),
        |_, pool, _| {
            let mut pool = pool.iter().copied().collect::<Vec<_>>();
            pool.sort_unstable();
            Ok(result(
                pool.chunks(2).map(<[usize]>::to_vec).collect(),
                Vec::new(),
            ))
        },
    );

    assert_eq!(ledger.promoted, 0);
    assert_eq!(ledger.refused_small, 5);
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
        &[],
        |pool, _, _| empty_knn(pool),
        |_, _, _| Ok(result(vec![vec![10, 11, 12, 13, 14, 15]], Vec::new())),
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
        &[],
        |pool, _, _| empty_knn(pool),
        |_, _, _| Ok(result(vec![vec![6, 7, 10, 11, 14, 15]], Vec::new())),
    );

    assert_eq!(
        map[&0],
        vec![0, 1, 2, 3, 4, 5],
        "a bin no proposal touched comes back whole"
    );
    assert_eq!(map[&1], vec![8, 9]);
    assert_eq!(map[&2], vec![12, 13]);
    assert!(map.values().any(|bin| bin == &vec![6, 7, 10, 11, 14, 15]));
    assert!(unbinned.is_empty());

    assert_eq!(ledger.pool_contigs, 16);
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
        bars: Bars {
            duplication_bar: 0.05,
            ..settings().bars
        },
        ..settings()
    };

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        settings,
        &[],
        |pool, _, _| empty_knn(pool),
        |_, _, _| Ok(result(vec![vec![0, 1, 2, 3], vec![4]], Vec::new())),
    );

    assert_eq!(ledger.dissolved_duplicated, 1, "{ledger}");
    assert_eq!(ledger.dissolved_small, 0);
    assert_eq!(ledger.refused_duplicated, 1);
    assert_eq!(ledger.promoted, 0);
    assert_eq!(map, start);
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
        settings(),
        &[],
        |pool, _, _| empty_knn(pool),
        |_, _, _| Ok(result(vec![vec![10, 11, 12, 13, 14, 15]], Vec::new())),
    );

    assert_eq!(ledger.dissolved_clean, 1, "{ledger}");
    assert_eq!(ledger.pool_contigs, 16);
    assert_eq!(map[&0], clean, "the clean bin came back whole");
    assert_eq!(map[&1], vec![8, 9]);
}

/// A round that accepts nothing at the genome floor is where the loop would stop, so the bar
/// drops a rung at a time and takes the same clusters rather than leaving them in the pool.
#[test]
fn the_ladder_takes_clusters_the_fixed_floor_refuses() {
    let (coverage, tnf, lengths) = pieces(16);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let clusters = vec![(4..9).collect::<Vec<_>>(), (9..14).collect::<Vec<_>>()];
    let mut map = BTreeMap::from([(0usize, vec![0, 1]), (1usize, vec![2, 3])]);
    let mut unbinned = (4..16).collect::<Vec<_>>();

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        settings(),
        &[],
        |pool, _, _| empty_knn(pool),
        |_, _, _| Ok(result(clusters.clone(), Vec::new())),
    );

    assert_eq!(
        ledger.promoted, 2,
        "500 kb clusters under a 600 kb genome floor: {ledger}"
    );
    assert!(ledger.rung > 0);
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
        rounds: 2,
        ..settings()
    };

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        settings,
        &[],
        |pool, _, _| empty_knn(pool),
        |_, _, round| {
            Ok(match round.n_neighbours {
                100 => result(vec![(0..12).collect()], Vec::new()),
                _ => result(vec![(0..6).collect(), (6..12).collect()], Vec::new()),
            })
        },
    );

    assert_eq!(ledger.rounds, 2 * POOL_VIEWS.len());
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
        rounds: 2,
        ..settings()
    };

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        settings,
        &[],
        |pool, _, _| empty_knn(pool),
        |_, _, round| {
            Ok(match round.n_neighbours {
                100 => result(vec![(0..8).collect()], Vec::new()),
                _ => result(vec![(8..16).collect()], Vec::new()),
            })
        },
    );

    assert_eq!(ledger.promoted, 2);
    assert_eq!(ledger.adopted_contigs, 16);
}

fn rungs(labellings: Vec<Vec<Vec<usize>>>) -> Vec<Partitioning> {
    labellings
        .into_iter()
        .map(|clusters| partitioning(clusters, Vec::new()))
        .collect()
}

/// The partition fits a ladder of labellings and its own score picks one. Reading only that one
/// hides every genome the winning rung merged, which is what the pool has a better judge for.
#[test]
fn every_labelling_handed_back_is_a_candidate() {
    let (coverage, tnf, lengths) = pieces(16);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let mut map = BTreeMap::from([(0usize, vec![0, 1])]);
    let mut unbinned = (2..16).collect::<Vec<_>>();

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        settings(),
        &[],
        |pool, _, _| empty_knn(pool),
        |_, _, _| Ok(rungs(vec![vec![(0..4).collect()], vec![(4..10).collect()]])),
    );

    assert_eq!(ledger.proposed, 2, "{ledger}");
    assert_eq!(
        ledger.promoted, 1,
        "the second rung holds the only cluster over the floor"
    );
    assert!(map.values().any(|bin| bin == &(4..10).collect::<Vec<_>>()));
}

/// One ranking over a fixed pool cannot see a genome the bigger one buries, because the graph
/// that buried it is never rebuilt.
#[test]
fn a_second_pass_searches_what_the_first_claimed_away() {
    let (coverage, tnf, lengths) = pieces(16);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let mut map = BTreeMap::from([(0usize, vec![0, 1])]);
    let mut unbinned = (2..16).collect::<Vec<_>>();
    let seen = std::cell::RefCell::new(Vec::new());

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        DissolveSettings {
            passes: 3,
            ..settings()
        },
        &[],
        |pool, _, _| empty_knn(pool),
        |_, pool, _| {
            let mut pool = pool.iter().copied().collect::<Vec<_>>();
            pool.sort_unstable();
            seen.borrow_mut().push(pool.len());
            Ok(result(
                vec![pool.iter().copied().take(6).collect()],
                Vec::new(),
            ))
        },
    );

    assert_eq!(ledger.promoted, 3, "{ledger}");
    let mut handed = seen.into_inner();
    handed.dedup();
    assert_eq!(
        handed,
        vec![16, 10, 4],
        "each pass is handed what the pass before it left"
    );
}

/// A pass finds worse bins than the one before it long before it finds none, so the budget is
/// a cap and the model's own scores say where to stop under it.
#[test]
fn the_passes_stop_once_a_pass_finds_worse_bins() {
    let (coverage, tnf, lengths) = pieces(32);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let mut map = BTreeMap::from([(0usize, vec![0, 1])]);
    let mut unbinned = (2..32).collect::<Vec<_>>();
    let seen = std::cell::RefCell::new(0usize);

    let ledger = dissolve(
        &features,
        None,
        &mut map,
        &mut unbinned,
        DissolveSettings {
            passes: 4,
            ..settings()
        },
        &[],
        |pool, _, _| empty_knn(pool),
        |_, pool, _| {
            let mut pool = pool.iter().copied().collect::<Vec<_>>();
            pool.sort_unstable();
            let pass = *seen.borrow();
            *seen.borrow_mut() += 1;
            let take = 8usize.saturating_sub(pass * 2).max(6);
            Ok(result(
                vec![pool.iter().copied().take(take).collect()],
                Vec::new(),
            ))
        },
    );

    assert_eq!(
        ledger.passes, 2,
        "the second pass is worth less than the first, so there is no third: {ledger}"
    );
}

/// The probe only reads if the group actually reaches the heap, so a group the search never
/// proposes has to be adoptable on its own.
#[test]
fn an_oracle_group_the_search_never_proposes_is_still_taken() {
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
        &[vec![0, 3, 6, 9]],
        |pool, _, _| empty_knn(pool),
        |_, _, _| Ok(result(Vec::new(), Vec::new())),
    );

    assert_eq!(ledger.promoted, 1);
    assert!(map.values().any(|bin| bin == &vec![0, 3, 6, 9]), "{map:?}");
}
