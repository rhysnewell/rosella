use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::learned::optim::{Plateau, Schedule};

fn stops_after(fall: f64, patience: usize) -> Option<usize> {
    let mut rng = StdRng::seed_from_u64(7);
    let mut plateau = Plateau::new(patience);
    (0..40).find(|epoch| {
        for _ in 0..50 {
            plateau.record(5.0 - fall * *epoch as f64 + rng.random_range(-0.1..0.1));
        }
        plateau.close().1
    })
}

#[test]
fn a_flat_loss_stops_and_a_falling_one_trains_on() {
    assert!(stops_after(0.0, 2).is_some_and(|epoch| epoch <= 3));
    assert!(stops_after(0.0, 5).is_some_and(|epoch| (5..=6).contains(&epoch)));
    assert_eq!(stops_after(0.0, 0), None);
    assert_eq!(stops_after(0.05, 2), None);
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
