//! Putting the bins under the floor back in the pot and embedding them again.

use std::collections::{BTreeMap, HashMap, HashSet};

use ndarray::Array2;
use rosella::clustering::clusterer::HDBSCANResult;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::rescue::{RescueSettings, rescue};

#[path = "../support/sketches.rs"]
mod sketches;

use sketches::{Fixture, grow, sibling};

const PIECE: usize = 100_000;
const FLOOR: usize = 200_000;
const GENOME: usize = 600_000;

fn settings() -> RescueSettings {
    RescueSettings {
        min_bin_size: FLOOR,
        genome_floor: Some(GENOME),
        duplication_bar: 1.0,
        min_contigs: 3,
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

    let ledger = rescue(&features, &mut map, &mut unbinned, settings(), |pool| {
        Ok(result(
            pool.iter().map(|contig| vec![*contig]).collect(),
            Vec::new(),
        ))
    });

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

    let ledger = rescue(&features, &mut map, &mut unbinned, settings(), |_| {
        Ok(result(vec![vec![10, 11, 12, 13, 14, 15]], Vec::new()))
    });

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

    let ledger = rescue(&features, &mut map, &mut unbinned, settings(), |_| {
        Ok(result(vec![vec![6, 7, 10, 11, 14, 15]], Vec::new()))
    });

    assert_eq!(map[&0], vec![0, 1, 2, 3, 4, 5], "a whole bin never dissolves");
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

    let ledger = rescue(&features, &mut map, &mut unbinned, settings(), |_| {
        Ok(result(
            vec![vec![6, 7, 8, 9, 10, 11], vec![10, 11, 12, 13, 14, 15]],
            Vec::new(),
        ))
    });

    assert_eq!(ledger.proposed, 2);
    assert_eq!(ledger.promoted, 0);
    assert_eq!(map, start);
    assert_eq!(unbinned, vec![10, 11, 12, 13, 14, 15]);
}

fn fused() -> Fixture {
    let native = grow(3, 300_000);
    Fixture::new(
        "rescue_fused",
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
    let settings = RescueSettings {
        min_bin_size: FLOOR,
        genome_floor: Some(300_000),
        duplication_bar: 0.05,
        min_contigs: 3,
    };

    let ledger = rescue(&features, &mut map, &mut unbinned, settings, |_| {
        Ok(result(vec![vec![0, 1, 2, 3], vec![4]], Vec::new()))
    });

    assert_eq!(ledger.dissolved_duplicated, 1, "{ledger}");
    assert_eq!(ledger.dissolved_small, 0);
    assert_eq!(ledger.refused_duplicated, 1);
    assert_eq!(ledger.refused_small, 1);
    assert_eq!(ledger.promoted, 0);
    assert_eq!(map, start);
}
