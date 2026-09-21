//! Whether a claim carved out of one bin shows the duplication that made the bin fail.

use std::collections::HashSet;

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::dissolve::Pot;

#[path = "../support/scorer.rs"]
mod scorer;

use scorer::GenomeScorer;

const PIECE: usize = 100_000;
const WORTH: f64 = 2.0;
const FLOOR: usize = 200_000;

fn features(lengths: &[usize]) -> (Array2<f64>, Array2<f64>) {
    (
        Array2::zeros((lengths.len(), 2)),
        Array2::zeros((lengths.len(), 2)),
    )
}

/// One genome across ten contigs: every piece of it reads clean, so nothing in the worth or
/// the bars refuses the half, and the parent has no duplication for the half to explain.
#[test]
fn half_of_one_organism_explains_nothing() {
    let lengths = vec![PIECE; 10];
    let (coverage, tnf) = features(&lengths);
    let held = ContigFeatures::new(&coverage, &tnf, &lengths);
    let quality = GenomeScorer::new(vec![0; 10]);

    let dissolved = vec![(0usize, (0..10).collect::<Vec<_>>())];
    let pot = Pot::new(&held, &quality, WORTH, FLOOR, &dissolved);
    let pool = (0..10).collect::<HashSet<_>>();

    let half = (0..5).collect::<Vec<_>>();
    assert!(
        !pot.unifies(&half, &pool, &HashSet::new()),
        "the bin carries one copy of the markers and the half leaves five contigs standing"
    );
    assert!(
        pot.unifies(&half, &pool, &(5..10).collect()),
        "nothing stands once the rest of the bin is claimed"
    );
}

/// Two genomes in one bin, which is what the pool is for. Taking one of them leaves the
/// other, and the claim reads clean where the bin it came from read doubled.
#[test]
fn an_organism_pulled_out_of_a_doubled_bin_is_taken() {
    let lengths = vec![PIECE; 10];
    let (coverage, tnf) = features(&lengths);
    let held = ContigFeatures::new(&coverage, &tnf, &lengths);
    let quality = GenomeScorer::new(vec![0, 0, 0, 0, 0, 1, 1, 1, 1, 1]);

    let dissolved = vec![(0usize, (0..10).collect::<Vec<_>>())];
    let pot = Pot::new(&held, &quality, WORTH, FLOOR, &dissolved);
    let pool = (0..10).collect::<HashSet<_>>();

    assert!(
        pot.unifies(&(0..5).collect::<Vec<_>>(), &pool, &HashSet::new()),
        "one genome of the two leaves the other behind"
    );
}

/// A claim drawing on several bins is a merge whatever it leaves behind, so the carve test
/// never reaches it.
#[test]
fn a_claim_drawing_on_two_bins_is_never_a_carve() {
    let lengths = vec![PIECE; 12];
    let (coverage, tnf) = features(&lengths);
    let held = ContigFeatures::new(&coverage, &tnf, &lengths);
    let quality = GenomeScorer::new(vec![0; 12]);

    let dissolved = vec![
        (0usize, (0..6).collect::<Vec<_>>()),
        (1usize, (6..12).collect::<Vec<_>>()),
    ];
    let pot = Pot::new(&held, &quality, WORTH, FLOOR, &dissolved);
    let pool = (0..12).collect::<HashSet<_>>();

    assert!(
        pot.unifies(&vec![0, 6], &pool, &HashSet::new()),
        "one contig from each bin is a join"
    );
}
