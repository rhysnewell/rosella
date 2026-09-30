use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::fillers::eject_fillers;

use crate::scorer::MarkerScorer;

const PIECE: usize = 100_000;

// Contigs 0 to 3 are the parent's, 4 to 6 come from another bin. 4 and 5 bring markers the
// parent lacks and 6 brings none.
fn markers() -> Vec<Vec<usize>> {
    vec![
        (0..4).collect(),
        (4..8).collect(),
        (8..12).collect(),
        (12..16).collect(),
        vec![16],
        vec![17],
        vec![],
    ]
}

fn eject(depths: &[f64]) -> (Vec<usize>, Vec<Vec<usize>>) {
    let mut coverage = Array2::zeros((depths.len(), 2));
    for (contig, depth) in depths.iter().enumerate() {
        coverage[[contig, 0]] = *depth;
    }
    let tnf = Array2::zeros((depths.len(), 2));
    let lengths = vec![PIECE; depths.len()];
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let quality = MarkerScorer::new(markers(), 20);
    let dissolved = vec![(0, vec![0, 1, 2, 3]), (1, vec![4, 5, 6])];
    let mut promoted = vec![vec![0, 1, 2, 3, 4, 5, 6]];
    let ejected = eject_fillers(&features, &quality, &dissolved, &mut promoted);
    (ejected, promoted)
}

#[test]
fn only_a_filler_off_the_parents_depth_goes() {
    let (ejected, promoted) = eject(&[9.0, 10.0, 11.0, 10.0, 2.0, 10.5, 2.0]);

    assert_eq!(ejected, vec![4]);
    assert_eq!(promoted, vec![vec![0, 1, 2, 3, 5, 6]]);
}

#[test]
fn a_parent_holding_an_uncovered_contig_refuses_nothing() {
    let (ejected, _) = eject(&[0.0, 10.0, 11.0, 10.0, 2.0, 10.5, 2.0]);

    assert!(ejected.is_empty());
}
