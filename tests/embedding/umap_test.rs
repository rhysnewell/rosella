//! flight derives UMAP's curve from assembly contiguity rather than exposing it, so the
//! derivation is pinned against flight 1.7.0 rather than left to drift.

use rosella::embedding::umap::{curve_params, n_components, n_x};

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
fn curve_matches_flight() {
    let curve = curve_params(&LENGTHS);
    assert!((curve.a - 1.5414973).abs() < 1e-6, "a was {}", curve.a);
    assert!((curve.b - 0.3).abs() < 1e-6, "b was {}", curve.b);
}

#[test]
fn curve_stays_inside_flights_clamps() {
    for lengths in [vec![1500; 50], vec![50_000_000; 3], vec![1500, 9_000_000]] {
        let curve = curve_params(&lengths);
        assert!((1.4..=2.0).contains(&curve.a), "a was {}", curve.a);
        assert!((0.3..=0.4).contains(&curve.b), "b was {}", curve.b);
    }
}

#[test]
fn n_x_handles_an_empty_assembly() {
    assert_eq!(n_x(&[], 50.0), 0);
}

#[test]
fn components_track_sample_count_within_flights_bounds() {
    assert_eq!(n_components(1), 2);
    assert_eq!(n_components(4), 4);
    assert_eq!(n_components(50), 10);
}
