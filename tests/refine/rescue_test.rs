//! Putting the bins under the floor back in the pot and embedding them again.

use std::collections::{BTreeMap, HashMap, HashSet};

use ndarray::Array2;
use rosella::clustering::clusterer::HDBSCANResult;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::rescue::{RescueSettings, rescue};

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

/// Eight 100 kb contigs: one bin already over the genome floor, two under it, two loose.
fn fixture() -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    (
        Array2::zeros((10, 2)),
        Array2::zeros((10, 2)),
        vec![PIECE; 10],
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
    let (coverage, tnf, lengths) = fixture();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let mut map = bins();
    let mut unbinned = vec![9];

    let promoted = rescue(&features, &mut map, &mut unbinned, settings(), |pool| {
        Ok(result(
            pool.iter().map(|contig| vec![*contig]).collect(),
            Vec::new(),
        ))
    });

    assert_eq!(promoted, 0);
    assert_eq!(map, bins());
    assert_eq!(unbinned, vec![9]);
}

/// A bin dissolved into the pool and then not claimed by any accepted cluster has to come back
/// whole. Dropping it strands its long contigs as singletons for nothing.
#[test]
fn a_dissolved_bin_the_pool_does_not_claim_comes_back() {
    let lengths = vec![PIECE; 16];
    let coverage = Array2::zeros((16, 2));
    let tnf = Array2::zeros((16, 2));
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let mut map = BTreeMap::from([
        (0usize, vec![0, 1, 2, 3, 4, 5]),
        (1usize, vec![6, 7]),
        (2usize, vec![8, 9]),
    ]);
    let mut unbinned = vec![10, 11, 12, 13, 14, 15];

    let promoted = rescue(&features, &mut map, &mut unbinned, settings(), |_| {
        Ok(result(vec![vec![10, 11, 12, 13, 14, 15]], Vec::new()))
    });

    assert_eq!(promoted, 1);
    assert_eq!(map[&1], vec![6, 7], "an unclaimed bin is restored");
    assert_eq!(map[&2], vec![8, 9]);
    assert!(
        map.values().any(|bin| bin == &vec![10, 11, 12, 13, 14, 15]),
        "the claimed contigs became a bin: {map:?}"
    );
    assert!(unbinned.is_empty());
}
