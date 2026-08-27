//! Golden values generated from flight 1.7.0's numba metrics, so a divergence in the
//! Rust port shows up as a test failure rather than a benchmark regression.

use rosella::embedding::metrics::{euclidean, metabat, rho};

const TOLERANCE: f64 = 1e-9;

/// Interleaved per-sample coverage mean and variance, three samples.
const COVERAGE: [[f64; 6]; 4] = [
        [4.0, 2.0, 10.0, 5.0, 0.5, 1.0],
        [4.2, 2.1, 9.5, 4.8, 0.6, 1.2],
        [50.0, 9.0, 1.0, 0.5, 20.0, 3.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
];

/// Stand-in for centre log ratio transformed tetranucleotide vectors.
const TNF: [[f64; 5]; 4] = [
        [0.1, -0.2, 0.3, -0.4, 0.05],
        [0.12, -0.18, 0.29, -0.41, 0.06],
        [-0.5, 0.6, -0.1, 0.2, -0.2],
        [1.3, -0.7, 0.9, -1.1, 0.4],
];

const FLIGHT_METABAT: [f64; 16] = [
        1.0, 0.06570923836094418, 0.9983332221055775, 0.5627088505820356,
        0.06570923836094418, 1.0, 0.9976135237470866, 0.5927321719265695,
        0.9983332221055775, 0.9976135237470866, 1.0, 0.7261688046452667,
        0.5627088505820356, 0.5927321719265695, 0.7261688046452667, 1.0,
];

const FLIGHT_RHO: [f64; 16] = [
        0.0, 0.0015391822257913937, 1.5811623246492987, 0.5479028697571744,
        0.0015391822257913937, 0.0, 1.5821630056415796, 0.5422488591528162,
        1.5811623246492987, 1.5821630056415796, 0.0, 1.592051905920519,
        0.5479028697571744, 0.5422488591528162, 1.592051905920519, 0.0,
];

const FLIGHT_EUCLIDEAN: [f64; 16] = [
        0.0, 0.033166247903553984, 1.2579745625409124, 1.631716887208072,
        0.033166247903553984, 0.0, 1.2588089608832629, 1.6206788701035133,
        1.2579745625409124, 1.2588089608832629, 0.0, 2.824889378365107,
        1.631716887208072, 1.6206788701035133, 2.824889378365107, 0.0,
];

fn assert_matches_flight(name: &str, rows: &[&[f64]], expected: &[f64], metric: fn(&[f64], &[f64]) -> f64) {
    for (i, a) in rows.iter().enumerate() {
        for (j, b) in rows.iter().enumerate() {
            let want = expected[i * rows.len() + j];
            let got = metric(a, b);
            assert!(
                (got - want).abs() < TOLERANCE,
                "{}({}, {}): flight gives {}, rosella gives {}",
                name, i, j, want, got
            );
        }
    }
}

#[test]
fn metabat_matches_flight() {
    let rows: Vec<&[f64]> = COVERAGE.iter().map(|row| row.as_slice()).collect();
    assert_matches_flight("metabat", &rows, &FLIGHT_METABAT, metabat);
}

#[test]
fn rho_matches_flight() {
    let rows: Vec<&[f64]> = TNF.iter().map(|row| row.as_slice()).collect();
    assert_matches_flight("rho", &rows, &FLIGHT_RHO, rho);
}

#[test]
fn euclidean_matches_flight() {
    let rows: Vec<&[f64]> = TNF.iter().map(|row| row.as_slice()).collect();
    assert_matches_flight("euclidean", &rows, &FLIGHT_EUCLIDEAN, euclidean);
}

/// flight divides by zero on two constant vectors. Rosella calls them identical instead.
#[test]
fn rho_survives_constant_vectors() {
    assert_eq!(rho(&[0.0; 5], &[0.0; 5]), 0.0);
    assert_eq!(rho(&[3.0; 5], &[3.0; 5]), 0.0);
}

/// flight skips samples where the two means are equal, so a row against itself has no
/// samples left to average and comes back maximally distant. Pinned because it is
/// surprising, not because it is desirable.
#[test]
fn metabat_self_distance_is_maximal() {
    assert_eq!(metabat(&COVERAGE[0], &COVERAGE[0]), 1.0);
}
