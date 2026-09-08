//! The arm that takes a contig out because the bin already holds its sequence.

use std::collections::BTreeMap;

use rosella::refine::duplication::{DuplicationSettings, eject_duplicated, intruders};

#[path = "../support/sketches.rs"]
mod sketches;

use sketches::{Fixture, grow, sibling};

fn settings() -> DuplicationSettings {
    DuplicationSettings {
        bar: 0.05,
        link: 0.5,
        min_hashes: 20,
    }
}

fn fused() -> Fixture {
    let native = grow(3, 300_000);
    Fixture::new(
        "fused",
        vec![
            native.clone(),
            sibling(&native[0..30_000]),
            sibling(&native[100_000..130_000]),
            sibling(&native[200_000..230_000]),
            grow(41, 50_000),
        ],
    )
}

#[test]
fn the_swallowed_side_leaves_and_the_host_stays() {
    let fixture = fused();
    let mut leaving = intruders(&fixture.features(), &[0, 1, 2, 3, 4], settings())
        .into_iter()
        .map(|(contig, _)| contig)
        .collect::<Vec<_>>();
    leaving.sort_unstable();
    assert_eq!(leaving, vec![1, 2, 3]);
}

#[test]
fn a_bin_of_unrelated_contigs_is_left_alone() {
    let fixture = Fixture::new(
        "unrelated",
        (0..5).map(|seed| grow(100 + seed, 60_000)).collect(),
    );
    assert!(intruders(&fixture.features(), &[0, 1, 2, 3, 4], settings()).is_empty());
}

/// Containment is mutual between two near complete siblings, so only mass says which is which.
#[test]
fn between_two_siblings_only_the_smaller_leaves() {
    let native = grow(17, 200_000);
    let fixture = Fixture::new("mutual", vec![native.clone(), sibling(&native[0..160_000])]);
    let leaving = intruders(&fixture.features(), &[0, 1], settings())
        .into_iter()
        .map(|(contig, _)| contig)
        .collect::<Vec<_>>();
    assert_eq!(leaving, vec![1]);
}

#[test]
fn a_bin_under_the_bar_is_never_examined() {
    let fixture = fused();
    let calm = DuplicationSettings {
        bar: 0.9,
        ..settings()
    };
    assert!(intruders(&fixture.features(), &[0, 1, 2, 3, 4], calm).is_empty());
}

#[test]
fn the_drain_stops_before_the_bin_falls_under_the_floor() {
    let fixture = fused();
    let mut bins = BTreeMap::from([(0usize, vec![0, 1, 2, 3, 4])]);
    let ejected = eject_duplicated(&fixture.features(), &mut bins, settings(), 400_000);
    assert_eq!(ejected.len(), 1, "only one 30 kb piece fits above the floor");
    assert_eq!(bins[&0].len(), 4);
}

