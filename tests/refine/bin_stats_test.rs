//! Golden values from flight's `metrics.get_averages` on the same fixture the metric tests
//! use, so a drift in the aggregation shows up here rather than as a benchmark regression.

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::embedding::metrics::{euclidean, metabat, rho};
use rosella::refine::bin_stats::{AGGREGATE, EUCLIDEAN, METABAT, RHO, Thresholds, bin_stats};

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

const FLIGHT_PER_CONTIG: [[f64; 4]; 4] = [
    [
        0.5422504370161857,
        0.7102014588774215,
        0.9742858992175129,
        0.6554183265064363,
    ],
    [
        0.5520183113448668,
        0.708650349006729,
        0.9708846929634434,
        0.6599270419388382,
    ],
    [
        0.907371850165977,
        1.585125745403799,
        1.7805576339297609,
        1.0187897507836727,
    ],
    [
        0.6272032757179572,
        0.8940678782768364,
        2.025761711892231,
        0.8163425547490729,
    ],
];

const FLIGHT_MEANS: [f64; 4] = [
    0.6572109685612466,
    0.9745113578911965,
    1.437872484500737,
    0.7876194184945051,
];

fn fixture() -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let coverage = Array2::from_shape_vec((4, 6), COVERAGE.concat()).unwrap();
    let tnf = Array2::from_shape_vec((4, 5), TNF.concat()).unwrap();
    (coverage, tnf, vec![100_000; 4])
}

#[test]
fn matches_flight_get_averages() {
    let (coverage, tnf, lengths) = fixture();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let stats = bin_stats(&features, &[0, 1, 2, 3], 42).unwrap();

    for (position, expected) in FLIGHT_PER_CONTIG.iter().enumerate() {
        for column in 0..4 {
            assert!(
                (stats.per_contig[position][column] - expected[column]).abs() < TOLERANCE,
                "contig {position} column {column}: {} != {}",
                stats.per_contig[position][column],
                expected[column]
            );
        }
    }

    for column in 0..4 {
        assert!((stats.mean[column] - FLIGHT_MEANS[column]).abs() < TOLERANCE);
    }
}

#[test]
fn a_bin_of_one_has_no_statistics() {
    let (coverage, tnf, lengths) = fixture();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    assert!(bin_stats(&features, &[0], 42).is_none());
}

/// Rows of two well separated groups, so the sampled estimate has real structure to
/// recover rather than a single blob where any estimate would look right.
fn wide_fixture(n: usize) -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let mut coverage = Array2::zeros((n, 2));
    let mut tnf = Array2::zeros((n, 4));
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };

    for row in 0..n {
        let group = (row % 2) as f64;
        coverage[[row, 0]] = 5.0 + group * 40.0 + next();
        coverage[[row, 1]] = 2.0 + next();
        for column in 0..4 {
            tnf[[row, column]] = group - 0.5 + next() * 0.2;
        }
    }

    (coverage, tnf, vec![50_000; n])
}

fn brute_force_means(features: &ContigFeatures, n: usize) -> [f64; 4] {
    let weight = features.weight();
    let mut totals = [0.0f64; 4];
    for i in 0..n {
        let mut row = [0.0f64; 4];
        for j in 0..n {
            if i == j {
                continue;
            }
            let md = metabat(features.coverage_row(i), features.coverage_row(j));
            let proportionality = rho(features.tnf_row(i), features.tnf_row(j));
            row[METABAT] += md;
            row[RHO] += proportionality;
            row[EUCLIDEAN] += euclidean(features.tnf_row(i), features.tnf_row(j));
            row[AGGREGATE] += (md.powf(weight) * proportionality.powf(1.0 - weight)).sqrt();
        }
        for column in 0..4 {
            totals[column] += row[column] / (n - 1) as f64;
        }
    }
    totals.map(|total| total / n as f64)
}

/// Past the exact limit every contig is scored against one shared sample. The estimate has
/// to stay close to the full answer and has to be the same on every run.
#[test]
fn the_sampled_path_tracks_the_exact_one() {
    let n = 2_001;
    let (coverage, tnf, lengths) = wide_fixture(n);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let indices = (0..n).collect::<Vec<_>>();

    let sampled = bin_stats(&features, &indices, 42).unwrap();
    let exact = brute_force_means(&features, n);

    for column in 0..4 {
        assert!(
            (sampled.mean[column] - exact[column]).abs() < 0.05,
            "column {column}: {} against {}",
            sampled.mean[column],
            exact[column]
        );
    }

    let repeat = bin_stats(&features, &indices, 42).unwrap();
    assert_eq!(sampled.per_contig, repeat.per_contig);
}

#[test]
fn thresholds_average_over_the_large_bins_only() {
    let (coverage, tnf, lengths) = fixture();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let stats = bin_stats(&features, &[0, 1, 2, 3], 42).unwrap();

    let ignored = Thresholds::from_bins(std::iter::once((999_999, &stats)));
    assert_eq!(ignored.mean, [0.0; 4]);

    let counted = Thresholds::from_bins(std::iter::once((2_000_000, &stats)));
    assert_eq!(counted.mean, stats.mean);
}
