use rosella::refine::rung::{Bars, RUNGS};

const MIN_BIN_SIZE: usize = 200_000;
const TOP: usize = 3_000_000;

fn bars() -> Bars {
    Bars {
        min_bin_size: MIN_BIN_SIZE,
        completeness: 90.0,
        contamination: 5.0,
        worth: rosella::quality::Worth { contamination: 5.0, allowance: 0.0 },
        rung_floor: 0.56,
        rung_ceiling: 3.0,
    }
}

#[test]
fn a_scorer_blind_to_size_is_still_held_to_the_genome_scale() {
    assert_eq!(bars().at(TOP, 0, true).floor, MIN_BIN_SIZE);
    assert_eq!(bars().at(TOP, 0, false).floor, TOP);

    let floors = (0..RUNGS)
        .map(|rung| bars().at(TOP, rung, false).floor)
        .collect::<Vec<_>>();
    assert!(
        floors.windows(2).all(|pair| pair[0] >= pair[1]),
        "{floors:?}"
    );
    assert!(floors[RUNGS - 1] > MIN_BIN_SIZE, "{floors:?}");
}
