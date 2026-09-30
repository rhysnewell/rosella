use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::restore::{Judge, restore};
use rosella::refine::rung::Rung;

use crate::scorer::MarkerScorer;

const PIECE: usize = 100_000;
const CATALOGUE: usize = 20;

// Contigs 0 to 3 are a host holding 16 of the 20 markers, 4 is the host with none, 5 and 6
// repeat eight of the host's markers from another genome, and 7 is a bare pool contig.
fn markers() -> Vec<Vec<usize>> {
    vec![
        (0..4).collect(),
        (4..8).collect(),
        (8..12).collect(),
        (12..16).collect(),
        vec![],
        (0..4).collect(),
        (4..8).collect(),
        vec![],
    ]
}

fn rung(completeness: f64, contamination: f64) -> Rung {
    Rung {
        floor: 1,
        completeness,
        contamination,
        ..Default::default()
    }
}

fn run(
    dissolved: &[(usize, Vec<usize>)],
    promoted: Vec<Vec<usize>>,
) -> rosella::refine::restore::Restored {
    let contigs = markers().len();
    let coverage = Array2::zeros((contigs, 2));
    let tnf = Array2::zeros((contigs, 2));
    let lengths = vec![PIECE; contigs];
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let quality = MarkerScorer::new(markers(), CATALOGUE);
    let held = Judge {
        features: &features,
        quality: &quality,
        reported: rung(50.0, f64::INFINITY),
        accept: rung(90.0, 5.0),
    };
    restore(&held, 2.0, dissolved, promoted)
}

#[test]
fn a_bin_handed_back_level_on_markers_goes_back() {
    let held = run(&[(0, vec![0, 1, 2, 3, 4])], vec![vec![0, 1, 2, 3, 7]]);

    assert_eq!(held.bins, 1);
    assert!(held.promoted.is_empty());
}

#[test]
fn a_cleanup_past_its_spread_stands() {
    let promoted = vec![vec![0, 1, 2, 3]];
    let held = run(&[(0, vec![0, 1, 2, 3, 5, 6])], promoted.clone());

    assert_eq!(held.bins, 0);
    assert_eq!(held.promoted, promoted);
}
