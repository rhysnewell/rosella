//! Whether the centroid cut takes a fused bin apart on the genome boundary and leaves a
//! pure one whole, and whether the same test judges pieces someone else proposed.

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::bisect::{candidate, separates};

const CLOUD: usize = 30;
const CONTIG_LENGTH: usize = 20_000;
const MIN_BIN_SIZE: usize = 200_000;
const ELIGIBLE: usize = 50;
const FIRST: [f64; 6] = [0.1, -0.2, 0.3, -0.4, 0.2, -0.1];
const SECOND: [f64; 6] = [-0.3, 0.4, -0.1, 0.2, -0.5, 0.3];

fn clouds(bases: &[[f64; 6]]) -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let n = CLOUD * bases.len();
    let mut coverage = Array2::zeros((n, 2));
    let mut tnf = Array2::zeros((n, 6));
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut noise = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };
    for (cloud, base) in bases.iter().enumerate() {
        for row in cloud * CLOUD..(cloud + 1) * CLOUD {
            coverage[[row, 0]] = 10.0 + noise();
            coverage[[row, 1]] = 2.0;
            for (column, value) in base.iter().enumerate() {
                tnf[[row, column]] = value + 0.1 * noise();
            }
        }
    }
    (coverage, tnf, vec![CONTIG_LENGTH; n])
}

#[test]
fn two_genomes_at_one_depth_come_apart_on_composition() {
    let (coverage, tnf, lengths) = clouds(&[FIRST, SECOND]);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let indices = (0..2 * CLOUD).collect::<Vec<_>>();

    let pieces = candidate(&features, &indices, MIN_BIN_SIZE, ELIGIBLE, 42).expect("a cut");
    let mut sides = pieces
        .iter()
        .map(|piece| piece.iter().map(|i| i / CLOUD).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    sides.sort();
    assert!(sides[0].iter().all(|cloud| *cloud == 0), "{pieces:?}");
    assert!(sides[1].iter().all(|cloud| *cloud == 1), "{pieces:?}");
    assert_eq!(sides[0].len() + sides[1].len(), 2 * CLOUD);
}

#[test]
fn one_genome_stays_whole() {
    let (coverage, tnf, lengths) = clouds(&[FIRST, FIRST]);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let indices = (0..2 * CLOUD).collect::<Vec<_>>();

    assert!(candidate(&features, &indices, MIN_BIN_SIZE, ELIGIBLE, 42).is_none());
}

#[test]
fn proposed_pieces_are_judged_by_the_same_modes() {
    let (coverage, tnf, lengths) = clouds(&[FIRST, SECOND]);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let genomes = vec![(0..CLOUD).collect::<Vec<_>>(), (CLOUD..2 * CLOUD).collect()];
    assert!(separates(&features, &genomes, ELIGIBLE, 42));

    let halved = vec![
        (0..CLOUD / 2).collect::<Vec<_>>(),
        (CLOUD / 2..CLOUD).collect(),
    ];
    assert!(!separates(&features, &halved, ELIGIBLE, 42));
}
