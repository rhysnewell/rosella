//! flight derives UMAP's curve from assembly contiguity rather than exposing it, so the
//! derivation is pinned rather than left to drift. `a` is still flight 1.7.0's value. `b` is
//! not: its window moved to where CAMI I low and medium both put the optimum.

use rosella::embedding::umap::{curve_params, length_weights, n_components, n_x};

const LENGTHS: [usize; 20] = [
    1500, 1800, 2400, 3100, 4200, 5600, 7000, 9500, 12000, 15000, 21000, 34000, 58000, 90000,
    150000, 260000, 480000, 750000, 1200000, 2100000,
];

#[test]
fn n_x_matches_flight() {
    for (percent, expected) in [
        (10.0, 260000),
        (25.0, 750000),
        (50.0, 1200000),
        (75.0, 2100000),
        (90.0, 2100000),
    ] {
        assert_eq!(n_x(&LENGTHS, percent), expected, "at {}%", percent);
    }
}

#[test]
fn curve_pins_a_to_flight_and_b_to_the_measured_window() {
    let curve = curve_params(&LENGTHS);
    assert!((curve.a - 1.5414973).abs() < 1e-6, "a was {}", curve.a);
    assert!((curve.b - 0.6).abs() < 1e-6, "b was {}", curve.b);
}

#[test]
fn curve_stays_inside_its_clamps() {
    for lengths in [vec![1500; 50], vec![50_000_000; 3], vec![1500, 9_000_000]] {
        let curve = curve_params(&lengths);
        assert!((1.4..=2.0).contains(&curve.a), "a was {}", curve.a);
        assert!((0.5..=0.6).contains(&curve.b), "b was {}", curve.b);
    }
}

#[test]
fn n_x_handles_an_empty_assembly() {
    assert_eq!(n_x(&[], 50.0), 0);
}

#[test]
fn components_follow_the_estimate_and_fall_back_to_the_sample_count() {
    assert_eq!(n_components(Some(7.20), 1), 7);
    assert_eq!(n_components(Some(12.05), 5), 10);
    assert_eq!(n_components(Some(0.4), 4), 4);
    assert_eq!(n_components(None, 1), 2);
    assert_eq!(n_components(None, 50), 10);
}

#[test]
fn length_weights_are_empty_at_power_zero() {
    assert!(length_weights(&[1000, 2000, 4000], 0.0).is_empty());
    assert!(length_weights(&[], 1.0).is_empty());
}

/// The normalisation is what keeps the graph's maximum weight, and with it the threshold
/// that drops the weakest edges, where it was.
#[test]
fn length_weights_hold_a_geometric_mean_of_one() {
    let weights = length_weights(&[2000, 3000, 8000, 150_000], 1.0);
    let log_mean = weights.iter().map(|w| (*w as f64).ln()).sum::<f64>() / weights.len() as f64;
    assert!(
        log_mean.exp() - 1.0 < 1e-5,
        "geometric mean was {}",
        log_mean.exp()
    );
}

/// Unbounded, a megabase contig is drawn thousands of times more often than a 1.5 kb one
/// and the short contigs never move.
#[test]
fn length_weights_are_bounded_either_side() {
    let weights = length_weights(&[500, 3000, 20_000_000], 1.0);
    assert!(
        weights.iter().all(|w| (0.2..=5.0).contains(w)),
        "{weights:?}"
    );
    assert!(weights[0] < weights[1] && weights[1] < weights[2]);
}
