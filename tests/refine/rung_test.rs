use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::rung::{Bars, RUNGS, Verdict, judge};

use crate::scorer::GenomeScorer;

const MIN_BIN_SIZE: usize = 200_000;
const TOP: usize = 3_000_000;

fn bars() -> Bars {
    Bars {
        min_bin_size: MIN_BIN_SIZE,
        completeness: 90.0,
        contamination: 5.0,
        worth: 5.0,
        rung_floor: 0.56,
    }
}

#[test]
fn the_size_floor_falls_with_every_rung_but_never_to_the_bin_floor() {
    assert_eq!(bars().at(TOP, 0).floor, TOP);

    let floors = (0..RUNGS)
        .map(|rung| bars().at(TOP, rung).floor)
        .collect::<Vec<_>>();
    assert!(
        floors.windows(2).all(|pair| pair[0] >= pair[1]),
        "{floors:?}"
    );
    assert!(floors[RUNGS - 1] > MIN_BIN_SIZE, "{floors:?}");
}

#[test]
fn a_whole_clean_genome_under_the_genome_scale_is_adopted_but_not_under_the_bin_floor() {
    let lengths = vec![300_000, 300_000, 150_000];
    let coverage = Array2::zeros((lengths.len(), 2));
    let tnf = Array2::zeros((lengths.len(), 2));
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let scorer = GenomeScorer::new(vec![0, 0, 1]);
    let rung = bars().at(TOP, 0);

    assert_eq!(judge(&features, &scorer, &[0, 1], rung), Verdict::Adopt);
    assert_eq!(judge(&features, &scorer, &[2], rung), Verdict::TooSmall);
}
