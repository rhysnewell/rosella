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
        ladder: Default::default(),
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
fn a_shorter_ladder_still_lands_its_last_rung_on_the_floor_share() {
    let mut bars = bars();
    bars.ladder.rungs = 3;
    let last = bars.at(TOP, bars.ladder.rungs - 1);

    assert!((last.completeness - bars.completeness * bars.rung_floor).abs() < 1e-9);
}

#[test]
fn a_capped_ladder_never_opens_past_the_tier_it_is_counted_at() {
    let open = bars();
    let mut capped = bars();
    capped.ladder.contamination_cap = 2.0;

    assert!(open.at(TOP, RUNGS - 1).contamination > open.tier());
    assert!(
        (0..RUNGS).all(|rung| capped.at(TOP, rung).contamination <= capped.tier()),
        "{:?}",
        (0..RUNGS)
            .map(|rung| capped.at(TOP, rung).contamination)
            .collect::<Vec<_>>()
    );
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
