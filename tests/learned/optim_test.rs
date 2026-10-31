use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::learned::optim::{Plateau, Schedule};

// A fixed offset per step stands for a batch's contig mix and the short last batch. Every epoch
// repeats it. The steps spread widely and the epoch means barely move.
fn stops_after(fall: f64, repeated: f64, warm: usize, patience: usize) -> Option<usize> {
    let mut rng = StdRng::seed_from_u64(7);
    let mut plateau = Plateau::new(patience);
    (0..40).find(|epoch| {
        for step in 0..50 {
            let offset = repeated
                * if step == 49 {
                    -2.0
                } else {
                    (step as f64).sin()
                };
            plateau.record(5.0 - fall * *epoch as f64 + offset + rng.random_range(-0.1..0.1));
        }
        plateau.close(*epoch >= warm).1
    })
}

#[test]
fn a_flat_loss_stops_and_a_falling_one_trains_on() {
    for (fall, repeated, warm, patience, stop) in [
        (0.0, 0.0, 0, 2, Some(2..=6)),
        (0.0, 0.0, 0, 5, Some(5..=10)),
        (0.0, 0.0, 12, 5, Some(17..=22)),
        (0.0, 0.0, 0, 0, None),
        (0.05, 0.0, 0, 2, None),
        (0.03, 0.3, 0, 5, None),
    ] {
        let at = stops_after(fall, repeated, warm, patience);
        assert_eq!(
            at.is_some(),
            stop.is_some(),
            "fall {fall} repeated {repeated} warm {warm} patience {patience} stopped at {at:?}"
        );
        if let (Some(at), Some(stop)) = (at, stop) {
            assert!(
                stop.contains(&at),
                "patience {patience} warm {warm} stopped at {at}"
            );
        }
    }
}

#[test]
fn an_early_stop_still_anneals_to_the_end_rate() {
    let mut cut = Schedule::new(1000, 2e-3);
    let full = Schedule::new(1000, 2e-3);
    cut.anneal_from(500, 100);
    assert_eq!(cut.at(499), full.at(499));
    let (rate, momentum) = cut.at(599);
    let (end, top) = full.at(999);
    assert!((rate - end).abs() < 1e-15, "{rate} against {end}");
    assert!((momentum - top).abs() < 1e-12);
    assert!(cut.at(550).0 < full.at(550).0);
}
