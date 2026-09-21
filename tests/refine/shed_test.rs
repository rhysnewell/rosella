use std::collections::BTreeMap;

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::shed::Split;

use rosella::markers::{ContigMarkers, Hit, MarkerSet};
use rosella::quality::Bars;
use rosella::refine::shed::shed;

const TABLE: &str = "model_name\tdomain\n\
                     alpha\tbac120\n\
                     beta\tbac120\n\
                     gamma\tbac120\n";

const SETS: &str = "set\tmedian_genome_bp\tmax_genome_bp\n\
                    bac\t300000\t0\n";

fn hit(marker: u16) -> Hit {
    Hit {
        marker,
        partial: false,
    }
}

fn markers(per_contig: Vec<Vec<Hit>>, lengths: Vec<usize>) -> ContigMarkers {
    ContigMarkers::new(per_contig, MarkerSet::parse(TABLE).with_scales(SETS)).with_lengths(lengths)
}

fn bin(contigs: &[usize]) -> BTreeMap<usize, Vec<usize>> {
    BTreeMap::from([(0, contigs.to_vec())])
}

fn loose() -> impl Fn(&[usize]) -> bool {
    |_: &[usize]| false
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

    assert_eq!(
        shed(
            &mut bins,
            &mut unbinned,
            &held,
            open(),
            0.0,
            &loose(),
            None
        ),
        1
    );
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

    assert_eq!(
        shed(
            &mut bins,
            &mut unbinned,
            &held,
            open(),
            0.0,
            &loose(),
            None
        ),
        0
    );
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

    assert_eq!(
        shed(
            &mut bins,
            &mut unbinned,
            &held,
            open(),
            0.0,
            &loose(),
            None
        ),
        1
    );
    assert_eq!(unbinned, vec![1]);
}

#[test]
fn a_bin_shed_empty_is_dropped() {
    let held = markers(vec![vec![hit(0)], vec![hit(0)]], vec![20_000, 20_000]);
    let mut bins = bin(&[0, 1]);
    let mut unbinned = Vec::new();

    shed(
        &mut bins,
        &mut unbinned,
        &held,
        open(),
        0.0,
        &loose(),
        None,
    );
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

    let traced = held.redundant_traced(&[0, 1, 2], 0.0);
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
    let held = ContigMarkers::new(per_contig, MarkerSet::parse(&table))
        .with_lengths([vec![100_000; 39], vec![20_000]].concat());
    let members = (0..40).collect::<Vec<_>>();

    let mut bins = bin(&members);
    let mut unbinned = Vec::new();
    let bars = Bars {
        completeness: 80.0,
        contamination: 5.0,
    };
    assert_eq!(
        shed(
            &mut bins,
            &mut unbinned,
            &held,
            bars,
            0.0,
            &loose(),
            None
        ),
        0
    );

    let mut bins = bin(&members);
    let mut unbinned = Vec::new();
    let bars = Bars {
        completeness: 80.0,
        contamination: 0.0,
    };
    assert_eq!(
        shed(
            &mut bins,
            &mut unbinned,
            &held,
            bars,
            0.0,
            &loose(),
            None
        ),
        1
    );
    assert_eq!(unbinned, vec![39]);
}

const CLOUD: usize = 30;
const CONTIG_LENGTH: usize = 20_000;
const MIN_BIN_SIZE: usize = 200_000;
const FIRST: [f64; 6] = [0.1, -0.2, 0.3, -0.4, 0.2, -0.1];
const SECOND: [f64; 6] = [-0.3, 0.4, -0.1, 0.2, -0.5, 0.3];

fn clouds(bases: &[[f64; 6]]) -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let n = CLOUD * bases.len();
    let mut coverage = Array2::zeros((n, 2));
    let mut tnf = Array2::zeros((n, 6));
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut noise = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };
    for (cloud, base) in bases.iter().enumerate() {
        for row in cloud * CLOUD..(cloud + 1) * CLOUD {
            coverage[[row, 0]] = 10.0 + noise();
            coverage[[row, 1]] = 2.0;
            for (column, value) in base.iter().enumerate() {
                tnf[[row, column]] = value + 0.1 * noise();
            }
        }
    }
    (coverage, tnf, vec![CONTIG_LENGTH; n])
}

/// One carrier in each cloud, so the walk nominates a victim in one and its twin in the other.
fn fused_markers(n: usize, lengths: &[usize]) -> ContigMarkers {
    let mut per_contig = vec![Vec::new(); n];
    per_contig[0] = vec![hit(0), hit(1)];
    per_contig[CLOUD] = vec![hit(0), hit(1)];
    per_contig[1] = vec![hit(2)];
    markers(per_contig, lengths.to_vec())
}

#[test]
fn a_bin_the_markers_call_fused_is_split_on_the_boundary_rather_than_thinned() {
    let (coverage, tnf, lengths) = clouds(&[FIRST, SECOND]);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let held = fused_markers(lengths.len(), &lengths);
    let mut bins = bin(&(0..2 * CLOUD).collect::<Vec<_>>());
    let mut unbinned = Vec::new();

    let split = Split {
        features: &features,
        min_bin_size: MIN_BIN_SIZE,
        seed: 42,
    };
    assert_eq!(
        shed(
            &mut bins,
            &mut unbinned,
            &held,
            open(),
            0.0,
            &loose(),
            Some(split)
        ),
        0
    );
    assert!(unbinned.is_empty());
    assert_eq!(bins.len(), 2);
    let mut sides = bins
        .values()
        .map(|piece| piece.iter().map(|index| index / CLOUD).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    sides.sort();
    assert!(sides[0].iter().all(|cloud| *cloud == 0), "{bins:?}");
    assert!(sides[1].iter().all(|cloud| *cloud == 1), "{bins:?}");
}

#[test]
fn one_cloud_falls_back_to_the_eviction() {
    let (coverage, tnf, lengths) = clouds(&[FIRST, FIRST]);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let held = fused_markers(lengths.len(), &lengths);
    let mut bins = bin(&(0..2 * CLOUD).collect::<Vec<_>>());
    let mut unbinned = Vec::new();

    let split = Split {
        features: &features,
        min_bin_size: MIN_BIN_SIZE,
        seed: 42,
    };
    assert_eq!(
        shed(
            &mut bins,
            &mut unbinned,
            &held,
            open(),
            0.0,
            &loose(),
            Some(split)
        ),
        1
    );
    assert_eq!(bins.len(), 1);
    assert_eq!(unbinned, vec![0]);
}

#[test]
fn a_contig_longer_than_the_bar_is_not_a_passenger() {
    let held = markers(
        vec![vec![hit(0), hit(1)], vec![hit(2)], vec![hit(0), hit(1)]],
        vec![900_000, 400_000, 250_000],
    );
    let mut bins = bin(&[0, 1, 2]);
    let mut unbinned = Vec::new();
    assert_eq!(shed(&mut bins, &mut unbinned, &held, open(), 2.0, &loose(), None), 0);
    assert_eq!(bins[&0], vec![0, 1, 2]);

    let mut bins = bin(&[0, 1, 2]);
    let mut unbinned = Vec::new();
    assert_eq!(shed(&mut bins, &mut unbinned, &held, open(), 3.0, &loose(), None), 1);
    assert_eq!(unbinned, vec![2]);
}

