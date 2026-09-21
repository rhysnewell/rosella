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
fn only_a_set_smaller_than_the_widest_moves_the_size_floor() {
    let rung = bars().at(TOP, 0);

    assert_eq!(rung.scaled_floor(1.0), rung.floor);
    assert_eq!(rung.scaled_floor(0.0), rung.floor);
    assert!(rung.scaled_floor(0.29) < rung.floor);
    assert!(rung.scaled_floor(0.29) >= MIN_BIN_SIZE);
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
