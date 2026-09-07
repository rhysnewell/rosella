//! Golden values generated from flight 1.7.0's numba metrics, so a divergence in the
//! Rust port shows up as a test failure rather than a benchmark regression.

use ndarray::Array2;
use rosella::embedding::bands::DepthBands;
use rosella::embedding::metrics::{
    AggregateMetric, Combination, CompositionMetric, CoverageAggregation, DistanceSettings,
    MIN_VAR, euclidean, metabat_with, prepared::PreparedAggregate, rho, variance_floor,
};

const COMPOSITION_METRICS: [CompositionMetric; 5] = [
    CompositionMetric::Rho,
    CompositionMetric::Cosine,
    CompositionMetric::Aitchison,
    CompositionMetric::Hellinger,
    CompositionMetric::TetraZ,
];

const TOLERANCE: f64 = 1e-9;

/// flight scored every sample, so the golden values only reproduce with the skip disabled.
const NO_SKIP: f64 = 0.0;

fn geometric(a: &[f64], b: &[f64]) -> f64 {
    metabat_with(
        a,
        b,
        MIN_VAR,
        MIN_VAR,
        CoverageAggregation::Geometric,
        NO_SKIP,
        None,
    )
    .0
}
const EPSILON: f64 = 1e-6;

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

/// Nothing to average is agreement, not distance, and it has to read that way whichever mode
/// would have done the averaging. flight returned the maximum here and split agreeing contigs.
#[test]
fn metabat_reads_an_empty_average_as_agreement() {
    let absent = [0.0, 0.0, 0.0, 0.0];
    for aggregation in [
        CoverageAggregation::Geometric,
        CoverageAggregation::Arithmetic,
        CoverageAggregation::Max,
    ] {
        assert_eq!(
            metabat_with(&[], &[], MIN_VAR, MIN_VAR, aggregation, NO_SKIP, None),
            (EPSILON, 0)
        );
        assert_eq!(
            metabat_with(&absent, &absent, MIN_VAR, MIN_VAR, aggregation, 0.01, None),
            (EPSILON, 0)
        );
    }
}

/// The skip has to fire on mutual absence and only on mutual absence, and the count it returns
/// is what reweights coverage against composition, so both halves are pinned together.
#[test]
fn mutual_absence_drops_a_sample_but_a_shallow_contig_keeps_its_own() {
    let deep = [50.0, 50.0, 0.0, 0.0, 40.0, 40.0];
    let shallow = [0.0, 0.0, 0.0, 0.0, 0.4, 0.4];
    let (_, scored) = metabat_with(
        &deep,
        &shallow,
        MIN_VAR,
        MIN_VAR,
        CoverageAggregation::Arithmetic,
        0.01,
        None,
    );
    assert_eq!(
        scored, 2,
        "the mutually absent sample is the only one to go"
    );

    // At 0.9 the deep contig is under its own bar in sample three, where the shallow one at 0.4
    // is over its. One bar shared across the pair would drop that sample.
    let (_, own_bars) = metabat_with(
        &deep,
        &shallow,
        MIN_VAR,
        MIN_VAR,
        CoverageAggregation::Arithmetic,
        0.9,
        None,
    );
    assert_eq!(
        own_bars, 2,
        "the shallow contig's presence is judged on its own scale"
    );
}

/// One sample in three agrees and the other two do not. The geometric mean calls the pair
/// close on the strength of the one, which is the behaviour the other two modes exist to
/// test against.
#[test]
fn aggregation_decides_how_much_one_agreeing_sample_is_worth() {
    let a = [4.0, 2.0, 10.0, 5.0, 0.5, 1.0];
    let b = [4.0, 2.0, 90.0, 5.0, 40.0, 1.0];
    let distance =
        |aggregation| metabat_with(&a, &b, MIN_VAR, MIN_VAR, aggregation, NO_SKIP, None).0;

    let geometric = distance(CoverageAggregation::Geometric);
    let arithmetic = distance(CoverageAggregation::Arithmetic);
    let max = distance(CoverageAggregation::Max);

    assert!(geometric < 0.05, "geometric was {geometric}");
    assert!(max > 0.9, "max was {max}");
    assert!(geometric < arithmetic && arithmetic < max);
}

#[test]
fn the_variance_floor_only_moves_when_asked_and_stays_bounded() {
    assert_eq!(variance_floor(1, 3000, false), MIN_VAR);
    assert_eq!(variance_floor(3000, 3000, true), MIN_VAR);
    assert_eq!(variance_floor(1, 3000, true), MIN_VAR * 2.0);
    assert_eq!(variance_floor(10_000_000, 3000, true), MIN_VAR * 0.25);
}

/// A sharper floor on a long contig has to make it more discriminating, not less.
#[test]
fn a_length_scaled_floor_separates_coverages_the_flat_floor_blurs() {
    let a = [4.0, 0.05, 10.0, 0.05];
    let b = [5.0, 0.05, 11.0, 0.05];
    let flat = metabat_with(
        &a,
        &b,
        MIN_VAR,
        MIN_VAR,
        CoverageAggregation::Geometric,
        NO_SKIP,
        None,
    )
    .0;
    let long = variance_floor(10_000_000, 3000, true);
    let sharp = metabat_with(
        &a,
        &b,
        long,
        long,
        CoverageAggregation::Geometric,
        NO_SKIP,
        None,
    )
    .0;
    assert!(sharp > flat, "flat {flat}, sharp {sharp}");
}

/// The prepared path stores the centred composition half as `f32`, so it agrees to single
/// precision rather than to the bit. Tight enough that a wrong term cannot hide under it.
#[test]
fn the_prepared_metric_agrees_with_the_pairwise_one() {
    let rows = (0..COVERAGE.len())
        .map(|row| {
            let mut values = COVERAGE[row].to_vec();
            values.extend_from_slice(&TNF[row]);
            values
        })
        .collect::<Vec<_>>();
    let floors = [MIN_VAR, 0.5, 2.0, 1.5];
    let coverage = Array2::from_shape_vec(
        (COVERAGE.len(), COVERAGE[0].len()),
        COVERAGE.iter().flatten().copied().collect(),
    )
    .expect("the fixture is rectangular");
    let drops = DepthBands::new(&coverage, 1, false);
    let keeps = DepthBands::new(&coverage, 1, true);

    for aggregation in [
        CoverageAggregation::Geometric,
        CoverageAggregation::Arithmetic,
        CoverageAggregation::Max,
    ] {
        for combination in [Combination::Geometric, Combination::Arithmetic] {
            for composition in COMPOSITION_METRICS {
                let settings = DistanceSettings {
                    aggregation,
                    combination,
                    composition,
                    composition_scale: 1.7,
                    presence_fraction: 0.01,
                    ..DistanceSettings::default()
                };
                for banded in [None, Some(&drops), Some(&keeps)] {
                    let pairwise =
                        AggregateMetric::new(COVERAGE[0].len(), settings).with_bands(banded);
                    let prepared =
                        PreparedAggregate::new(&rows, &floors, COVERAGE[0].len(), settings, banded);

                    for a in 0..rows.len() {
                        for b in 0..rows.len() {
                            let (left, right) = (
                                prepared.distance(a, b),
                                pairwise.distance(&rows[a], &rows[b], floors[a], floors[b]),
                            );
                            assert!(
                                (left - right).abs() <= 1e-5 * right.abs().max(1.0),
                                "{aggregation:?} {combination:?} {composition:?} disagreed on rows \
                         {a} and {b}: {left} against {right}"
                            );
                        }
                    }
                }
            }
        }
    }
}

/// The refiner's rho and aggregate bars are calibrated to [0, 2], so a variant that leaves the
/// range changes what every one of those thresholds means.
#[test]
fn every_composition_metric_stays_inside_the_calibrated_range() {
    let degenerate = vec![0.0; TNF[0].len()];
    for composition in COMPOSITION_METRICS {
        for a in 0..TNF.len() {
            let rows = TNF
                .iter()
                .map(|row| row.to_vec())
                .chain(std::iter::once(degenerate.clone()))
                .collect::<Vec<_>>();
            for b in 0..rows.len() {
                let distance = composition.distance(&rows[a], &rows[b], 1.7);
                assert!(
                    (0.0..=2.0).contains(&distance),
                    "{composition:?} left the range on rows {a} and {b}: {distance}"
                );
            }
        }
    }
}
