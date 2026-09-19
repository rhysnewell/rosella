use std::collections::BTreeMap;

use rosella::markers::{ContigMarkers, Hit, MarkerRules, MarkerSet};
use rosella::quality::Bars;
use rosella::refine::shed::shed;

const TABLE: &str = "model_name\tdomain\n\
                     alpha\tbac120\n\
                     beta\tbac120\n\
                     gamma\tbac120\n";

fn hit(marker: u16) -> Hit {
    Hit {
        marker,
        partial: false,
    }
}

fn markers(per_contig: Vec<Vec<Hit>>, lengths: Vec<usize>) -> ContigMarkers {
    ContigMarkers::new(per_contig, MarkerSet::parse(TABLE), MarkerRules::default())
        .with_lengths(lengths)
}

fn bin(contigs: &[usize]) -> BTreeMap<usize, Vec<usize>> {
    BTreeMap::from([(0, contigs.to_vec())])
}

fn open() -> Bars {
    Bars {
        completeness: f64::INFINITY,
        contamination: 0.0,
    }
}

#[test]
fn a_contig_whose_every_marker_the_bin_keeps_leaves() {
    let held = markers(
        vec![vec![hit(0), hit(1)], vec![hit(2)], vec![hit(0), hit(1)]],
        vec![900_000, 400_000, 20_000],
    );
    let mut bins = bin(&[0, 1, 2]);
    let mut unbinned = Vec::new();

    assert_eq!(shed(&mut bins, &mut unbinned, &held, open()), 1);
    assert_eq!(bins[&0], vec![0, 1]);
    assert_eq!(unbinned, vec![2]);
}

#[test]
fn the_last_carrier_of_a_marker_never_leaves() {
    let held = markers(
        vec![vec![hit(0), hit(0)], vec![hit(1)]],
        vec![10_000, 900_000],
    );
    let mut bins = bin(&[0, 1]);
    let mut unbinned = Vec::new();

    assert_eq!(shed(&mut bins, &mut unbinned, &held, open()), 0);
    assert!(unbinned.is_empty());
}

#[test]
fn shedding_one_copy_protects_the_other() {
    let held = markers(
        vec![vec![hit(0)], vec![hit(0)], vec![hit(1), hit(2)]],
        vec![30_000, 20_000, 900_000],
    );
    let mut bins = bin(&[0, 1, 2]);
    let mut unbinned = Vec::new();

    assert_eq!(shed(&mut bins, &mut unbinned, &held, open()), 1);
    assert_eq!(unbinned, vec![1]);
}

#[test]
fn a_bin_shed_empty_is_dropped() {
    let held = markers(vec![vec![hit(0)], vec![hit(0)]], vec![20_000, 20_000]);
    let mut bins = bin(&[0, 1]);
    let mut unbinned = Vec::new();

    shed(&mut bins, &mut unbinned, &held, open());
    assert_eq!(bins.len(), 1);
    assert_eq!(bins[&0].len(), 1);
}

#[test]
fn a_contig_carrying_a_marker_the_bin_lacks_completes_it() {
    let held = markers(
        vec![vec![hit(0)], vec![hit(1)], vec![hit(0)]],
        vec![900_000, 20_000, 20_000],
    );

    assert!(held.completes(&[0, 1], 1));
    assert!(!held.completes(&[0, 2], 2));
}

#[test]
fn the_trace_names_the_carrier_that_made_a_contig_look_redundant() {
    let held = markers(
        vec![vec![hit(0), hit(1)], vec![hit(2)], vec![hit(0), hit(1)]],
        vec![900_000, 400_000, 20_000],
    );

    let traced = held.redundant_traced(&[0, 1, 2]);
    assert_eq!(traced.len(), 1);
    assert_eq!(traced[0].contig, 2);
    assert_eq!(traced[0].markers, 2);
    assert_eq!(traced[0].twin, Some(0));
    assert_eq!(traced[0].shared, 2);
}

#[test]
fn a_bin_over_both_bars_keeps_its_duplicate() {
    let mut table = String::from("model_name\tdomain\n");
    for model in 0..40 {
        table.push_str(&format!("m{model}\tbac120\n"));
    }
    let mut per_contig = (0..39).map(|marker| vec![hit(marker)]).collect::<Vec<_>>();
    per_contig.push(vec![hit(0)]);
    let held = ContigMarkers::new(per_contig, MarkerSet::parse(&table), MarkerRules::default())
        .with_lengths([vec![100_000; 39], vec![20_000]].concat());
    let members = (0..40).collect::<Vec<_>>();

    let mut bins = bin(&members);
    let mut unbinned = Vec::new();
    let bars = Bars {
        completeness: 80.0,
        contamination: 5.0,
    };
    assert_eq!(shed(&mut bins, &mut unbinned, &held, bars), 0);

    let mut bins = bin(&members);
    let mut unbinned = Vec::new();
    let bars = Bars {
        completeness: 80.0,
        contamination: 0.0,
    };
    assert_eq!(shed(&mut bins, &mut unbinned, &held, bars), 1);
    assert_eq!(unbinned, vec![39]);
}
