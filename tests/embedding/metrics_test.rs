//! Golden values generated from flight 1.7.0's numba metrics, so a divergence in the
//! Rust port shows up as a test failure rather than a benchmark regression.

use rosella::embedding::metrics::{MIN_VAR, Moments, euclidean, metabat_with, overlap, rho};

const TOLERANCE: f64 = 1e-9;

/// flight scored every sample, so the golden values only reproduce with the skip disabled.
const NO_SKIP: f64 = 0.0;

const EPSILON: f64 = 1e-6;

/// The golden values are flight's geometric fold of the same per-sample overlaps this build
/// averages arithmetically, so the fold is done here and `overlap` stays the anchored half.
fn geometric(a: &[f64], b: &[f64]) -> f64 {
    let mut total = 0.0;
    let mut scored = 0usize;
    for (sample_a, sample_b) in a.chunks_exact(2).zip(b.chunks_exact(2)) {
        if sample_a[0] <= NO_SKIP && sample_b[0] <= NO_SKIP {
            continue;
        }
        let moments = |sample: &[f64]| Moments::new(sample[0], (sample[1] + EPSILON).max(MIN_VAR));
        total += overlap(moments(sample_a), moments(sample_b))
            .clamp(EPSILON, 1.0 - EPSILON)
            .ln();
        scored += 1;
    }
    match scored {
        0 => EPSILON,
        scored => (total / scored as f64).exp(),
    }
}

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

/// The diagonal is rosella's, not flight's: flight skipped agreeing samples and returned
/// 1.0 for a row against itself.
const FLIGHT_METABAT: [f64; 16] = [
    EPSILON,
    0.06570923836094418,
    0.9983332221055775,
    0.5627088505820356,
    0.06570923836094418,
    EPSILON,
    0.9976135237470866,
    0.5927321719265695,
    0.9983332221055775,
    0.9976135237470866,
    EPSILON,
    0.7261688046452667,
    0.5627088505820356,
    0.5927321719265695,
    0.7261688046452667,
    EPSILON,
];

const FLIGHT_RHO: [f64; 16] = [
    0.0,
    0.0015391822257913937,
    1.5811623246492987,
    0.5479028697571744,
    0.0015391822257913937,
    0.0,
    1.5821630056415796,
    0.5422488591528162,
    1.5811623246492987,
    1.5821630056415796,
    0.0,
    1.592051905920519,
    0.5479028697571744,
    0.5422488591528162,
    1.592051905920519,
    0.0,
];

const FLIGHT_EUCLIDEAN: [f64; 16] = [
    0.0,
    0.033166247903553984,
    1.2579745625409124,
    1.631716887208072,
    0.033166247903553984,
    0.0,
    1.2588089608832629,
    1.6206788701035133,
    1.2579745625409124,
    1.2588089608832629,
    0.0,
    2.824889378365107,
    1.631716887208072,
    1.6206788701035133,
    2.824889378365107,
    0.0,
];

fn assert_matches_flight(
    name: &str,
    rows: &[&[f64]],
    expected: &[f64],
    metric: fn(&[f64], &[f64]) -> f64,
) {
    for (i, a) in rows.iter().enumerate() {
        for (j, b) in rows.iter().enumerate() {
            let want = expected[i * rows.len() + j];
            let got = metric(a, b);
            assert!(
                (got - want).abs() < TOLERANCE,
                "{}({}, {}): flight gives {}, rosella gives {}",
                name,
                i,
                j,
                want,
                got
            );
        }
    }
}

#[test]
fn metabat_matches_flight() {
    let rows: Vec<&[f64]> = COVERAGE.iter().map(|row| row.as_slice()).collect();
    assert_matches_flight("metabat", &rows, &FLIGHT_METABAT, geometric);
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

/// Identical coverage is the strongest evidence two contigs share a genome, so it has to
/// come out closest. flight returned 1.0 here, its maximum.
#[test]
fn metabat_self_distance_is_minimal() {
    for row in COVERAGE.iter() {
        assert!(
            geometric(row, row) < 1e-5,
            "self distance was {}",
            geometric(row, row)
        );
    }
    assert!(geometric(&COVERAGE[0], &COVERAGE[0]) < geometric(&COVERAGE[0], &COVERAGE[1]));
}

/// Nothing to average is agreement, not distance. flight returned the maximum here and split
/// agreeing contigs.
#[test]
fn metabat_reads_an_empty_average_as_agreement() {
    let absent = [0.0, 0.0, 0.0, 0.0];
    assert_eq!(
        metabat_with(&[], &[], MIN_VAR, MIN_VAR, NO_SKIP),
        (EPSILON, 0)
    );
    assert_eq!(
        metabat_with(&absent, &absent, MIN_VAR, MIN_VAR, 0.01),
        (EPSILON, 0)
    );
}

/// The skip has to fire on mutual absence and only on mutual absence, and the count it returns
/// is what reweights coverage against composition, so both halves are pinned together.
#[test]
fn mutual_absence_drops_a_sample_but_a_shallow_contig_keeps_its_own() {
    let deep = [50.0, 50.0, 0.0, 0.0, 40.0, 40.0];
    let shallow = [0.0, 0.0, 0.0, 0.0, 0.4, 0.4];
    let (_, scored) = metabat_with(&deep, &shallow, MIN_VAR, MIN_VAR, 0.01);
    assert_eq!(
        scored, 2,
        "the mutually absent sample is the only one to go"
    );

    // At 0.9 the deep contig is under its own bar in sample three, where the shallow one at 0.4
    // is over its. One bar shared across the pair would drop that sample.
    let (_, own_bars) = metabat_with(&deep, &shallow, MIN_VAR, MIN_VAR, 0.9);
    assert_eq!(
        own_bars, 2,
        "the shallow contig's presence is judged on its own scale"
    );
}

/// One sample in three agrees and the other two do not. The arithmetic mean has to let the two
/// disagreeing samples carry the pair, which is what flight's geometric fold could not do.
#[test]
fn one_agreeing_sample_cannot_carry_a_disagreeing_pair() {
    let a = [4.0, 2.0, 10.0, 5.0, 0.5, 1.0];
    let b = [4.0, 2.0, 90.0, 5.0, 40.0, 1.0];
    let arithmetic = metabat_with(&a, &b, MIN_VAR, MIN_VAR, NO_SKIP).0;

    assert!(arithmetic > 0.6, "arithmetic was {arithmetic}");
    assert!(arithmetic > geometric(&a, &b));
}


/// The refiner's rho and aggregate bars are calibrated to [0, 2], so rho leaving that range
/// changes what every one of those thresholds means.
#[test]
fn rho_stays_inside_the_calibrated_range() {
    let rows = TNF
        .iter()
        .map(|row| row.to_vec())
        .chain(std::iter::once(vec![0.0; TNF[0].len()]))
        .collect::<Vec<_>>();
    for a in &rows {
        for b in &rows {
            let distance = rho(a, b);
            assert!((0.0..=2.0).contains(&distance), "rho left the range: {distance}");
        }
    }
}
