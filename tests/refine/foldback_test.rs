use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::foldback::fold_back;
use rosella::refine::rung::Rung;

use crate::scorer::MarkerScorer;

const PIECE: usize = 100_000;
const CATALOGUE: usize = 20;

// Contigs 0 to 3 hold 16 of the 20 markers, 4 fills marker 16, 5 repeats markers 0 to 3, 6 to 9
// are a second genome holding markers 0 to 11, and 10 carries none.
fn markers() -> Vec<Vec<usize>> {
    vec![
        (0..4).collect(),
        (4..8).collect(),
        (8..12).collect(),
        (12..16).collect(),
        vec![16],
        (0..4).collect(),
        (0..3).collect(),
        (3..6).collect(),
        (6..9).collect(),
        (9..12).collect(),
        vec![],
    ]
}

fn fold(dissolved: &[(usize, Vec<usize>)], mut promoted: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    let contigs = markers().len();
    let coverage = Array2::zeros((contigs, 2));
    let tnf = Array2::zeros((contigs, 2));
    let lengths = vec![PIECE; contigs];
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let quality = MarkerScorer::new(markers(), CATALOGUE);
    let reported = Rung {
        completeness: 50.0,
        contamination: f64::INFINITY,
        ..Default::default()
    };
    fold_back(&features, &quality, 2.0, reported, dissolved, &mut promoted);
    promoted
}

#[test]
fn only_the_shard_contig_that_fills_a_marker_comes_back() {
    let promoted = fold(&[(0, vec![0, 1, 2, 3, 4, 5, 10])], vec![vec![0, 1, 2, 3]]);

    assert_eq!(promoted, vec![vec![0, 1, 2, 3, 4]]);
}

#[test]
fn a_leftover_that_reports_as_a_bin_stays_apart() {
    let claim = vec![0, 1, 2, 3];
    let promoted = fold(&[(0, vec![0, 1, 2, 3, 4, 6, 7, 8, 9])], vec![claim.clone()]);

    assert_eq!(promoted, vec![claim]);
}

#[test]
fn a_claim_holding_half_the_bin_or_less_takes_nothing_back() {
    let claim = vec![0, 1];
    let promoted = fold(&[(0, vec![0, 1, 2, 3])], vec![claim.clone()]);

    assert_eq!(promoted, vec![claim]);
}
