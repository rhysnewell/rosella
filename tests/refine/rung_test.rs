use rosella::refine::rung::{Bars, RUNGS};

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
