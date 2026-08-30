use std::collections::BTreeMap;

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::merger::merge_bins;

const CONTIG_LENGTH: usize = 100_000;
const PER_BIN: usize = 6;

fn jitter(seed: u64) -> impl FnMut() -> f64 {
    let mut state = seed;
    move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn two_clouds() -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let rows = PER_BIN * 2;
    let mut coverage = Array2::zeros((rows, 2));
    let mut tnf = Array2::zeros((rows, 4));
    let mut next = jitter(0x2545_F491_4F6C_DD1D);

    for row in 0..rows {
        let second = if row < PER_BIN { 0.0 } else { 1.0 };
        coverage[[row, 0]] = 10.0 + second * 4.0 + 6.0 * next();
        coverage[[row, 1]] = 4.0 + next();
        for column in 0..4 {
            tnf[[row, column]] = 0.1 * column as f64 + second * 0.2 + 0.4 * next();
        }
    }
    (coverage, tnf, vec![CONTIG_LENGTH; rows])
}

/// Centroids 0.308 apart against a bar of 0.345, so the first gate accepts, and a union at
/// 0.395 that the second refuses. Removing that refusal costs CAMI I high 45 bins at t5, so
/// the rejection is load-bearing however badly the two quantities match.
#[test]
fn the_second_gate_refuses_what_the_first_accepts() {
    let (coverage, tnf, lengths) = two_clouds();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let (merged, merges) = merge_bins(
        &features,
        BTreeMap::from([
            (0, (0..PER_BIN).collect::<Vec<_>>()),
            (1, (PER_BIN..PER_BIN * 2).collect()),
        ]),
        usize::MAX,
        42,
    );

    assert_eq!(merges, 0);
    assert_eq!(merged.len(), 2);
}
