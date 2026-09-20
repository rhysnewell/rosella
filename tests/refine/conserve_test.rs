//! Whether a pool claim is allowed to leave the bins it drew from worse than it found them.

use std::collections::HashSet;

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::dissolve::Pot;

#[path = "../support/scorer.rs"]
mod scorer;

use scorer::GenomeScorer;

const PIECE: usize = 100_000;
const WORTH: f64 = 2.0;

#[test]
fn a_claim_that_leaves_a_bin_worse_is_refused_where_the_guard_takes_it() {
    let coverage = Array2::zeros((10, 2));
    let tnf = Array2::zeros((10, 2));
    let lengths = vec![PIECE; 10];
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let quality = GenomeScorer::new(vec![0; 10]);

    let dissolved = vec![
        (0usize, vec![0, 1, 2, 3, 4, 5, 6]),
        (1usize, vec![7, 8, 9]),
    ];
    let pot = Pot::new(&features, &quality, WORTH, &dissolved);
    let pool = (0..10).collect::<HashSet<_>>();
    let claimed = HashSet::new();

    let carve = vec![0, 1, 2, 7, 8, 9];
    assert!(
        pot.takeable(&carve),
        "the shipped guard takes a minority of bin 0 untested"
    );
    assert!(
        !pot.conserves(&carve, &pool, &claimed),
        "bin 0 falls from 7 contigs to 4 and neither piece is worth what it was"
    );

    let unite = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
    assert!(
        pot.conserves(&unite, &pool, &claimed),
        "a claim that takes both bins whole is worth more than either"
    );
}

/// What is left of a bin shrinks as a pass claims from it, and the test is against the bin as
/// it stands rather than as it was dissolved, or the first claim spends the whole allowance.
#[test]
fn the_bin_a_claim_is_weighed_against_is_the_one_still_standing() {
    let coverage = Array2::zeros((10, 2));
    let tnf = Array2::zeros((10, 2));
    let lengths = vec![PIECE; 10];
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let quality = GenomeScorer::new(vec![0; 10]);

    let dissolved = vec![(0usize, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9])];
    let pot = Pot::new(&features, &quality, WORTH, &dissolved);
    let pool = (0..10).collect::<HashSet<_>>();

    let claim = vec![0, 1, 2, 3, 4];
    assert!(
        !pot.conserves(&claim, &pool, &HashSet::new()),
        "half of a whole bin is worth less than the whole"
    );
    assert!(
        pot.conserves(&claim, &pool, &(5..10).collect()),
        "once the rest is already claimed the standing bin is what is left"
    );
}
