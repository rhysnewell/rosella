use rosella::embedding::metrics::calibration::{LengthCalibration, sample_pairs};

const SLOPE: f64 = 120.0;
const FLOOR: f64 = 0.05;

fn synthetic() -> Vec<(f64, f64)> {
    let mut pairs = Vec::new();
    let mut state = 1u64;
    let draw = |state: &mut u64| {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        (*state >> 40) as f64 / (1u64 << 24) as f64
    };
    for _ in 0..4_000 {
        let length = |unit: f64| 1_500.0 * (500_000.0f64 / 1_500.0).powf(unit);
        let sum = 1.0 / length(draw(&mut state)) + 1.0 / length(draw(&mut state));
        let spread = draw(&mut state) * 0.8;
        pairs.push((sum, SLOPE * sum + FLOOR + spread));
    }
    pairs
}

#[test]
fn the_fit_recovers_the_length_term_from_the_lower_edge() {
    let calibration = LengthCalibration::fit(&synthetic()).expect("enough pairs to fit");
    let short = 2.0 / 1_500.0;
    let long = 2.0 / 200_000.0;

    let raw_short = SLOPE * short + FLOOR;
    let raw_long = SLOPE * long + FLOOR;
    assert!(
        raw_short - raw_long > 0.1,
        "the untouched gap is worth closing"
    );

    let gap = calibration.apply(raw_short, short) - calibration.apply(raw_long, long);
    let closed = 1.0 - gap.abs() / (raw_short - raw_long);
    assert!(
        closed > 0.8,
        "the length term should mostly come out, closed {closed:.2} of it"
    );
}

#[test]
fn a_flat_relationship_leaves_the_distance_where_it_was() {
    let pairs = (0..500)
        .map(|step| (1.0 / (2_000.0 + step as f64 * 100.0), 0.4))
        .collect::<Vec<_>>();
    let calibration = LengthCalibration::fit(&pairs).expect("enough pairs to fit");
    let moved = calibration.apply(0.9, 1.0 / 5_000.0) - 0.9;
    assert!(moved.abs() < 1e-9, "nothing to correct, moved by {moved}");
}

#[test]
fn too_few_pairs_refuses_rather_than_fitting_noise() {
    assert!(LengthCalibration::fit(&[(0.001, 0.4), (0.002, 0.5)]).is_none());
}

#[test]
fn the_sample_is_the_same_on_every_call() {
    let mut first = Vec::new();
    let mut second = Vec::new();
    sample_pairs(400, |a, b| first.push((a, b)));
    sample_pairs(400, |a, b| second.push((a, b)));
    assert_eq!(first, second);
    assert!(first.iter().all(|(a, b)| a < b && *b < 400));
}
