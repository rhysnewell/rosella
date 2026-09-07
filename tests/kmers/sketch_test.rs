use std::collections::HashSet;

use rosella::kmers::sketch::{SketchParams, sketch_sequence};

const K: u8 = 31;

fn params(scale: u64) -> SketchParams {
    SketchParams {
        kmer_size: K,
        scale,
    }
}

fn grow(seed: u64, length: usize) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 33) as usize % 4]
        })
        .collect()
}

fn drift(sequence: &[u8], every: usize) -> Vec<u8> {
    sequence
        .iter()
        .enumerate()
        .map(|(at, base)| if at % every == 0 { b'A' + (base % 3) } else { *base })
        .collect()
}

fn duplication(runs: &[(Vec<u64>, u32)]) -> f64 {
    let total: u32 = runs.iter().map(|(_, occurrences)| occurrences).sum();
    let distinct = runs
        .iter()
        .flat_map(|(hashes, _)| hashes.iter().copied())
        .collect::<HashSet<_>>();
    1.0 - distinct.len() as f64 / f64::from(total)
}

/// The whole point of sketching is that the scale divides out of both halves of the
/// duplication fraction, so an exact count and a 200-fold subsample have to agree.
#[test]
fn the_duplication_fraction_survives_subsampling() {
    let native = grow(11, 40000);
    let sibling = drift(&native, 100);
    let exact = duplication(&[
        sketch_sequence(&native, params(1)),
        sketch_sequence(&sibling, params(1)),
    ]);
    let sketched = duplication(&[
        sketch_sequence(&native, params(200)),
        sketch_sequence(&sibling, params(200)),
    ]);
    assert!(exact > 0.2, "the fixture has to be visibly duplicated, got {exact}");
    assert!(
        (exact - sketched).abs() < 0.03,
        "exact {exact} against sketched {sketched}"
    );
}

/// A fold onto the wrong strand makes every bin look pure, which no scored run would flag.
#[test]
fn a_strand_and_its_complement_sketch_alike() {
    let sequence = grow(29, 20000);
    let flipped = sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            _ => b'C',
        })
        .collect::<Vec<u8>>();
    let (forward, _) = sketch_sequence(&sequence, params(50));
    let (reverse, _) = sketch_sequence(&flipped, params(50));
    assert_eq!(forward, reverse);
}

#[test]
fn an_ambiguous_base_does_not_join_the_flanks() {
    let left = grow(5, 4000);
    let right = grow(7, 4000);
    let mut spanning = left.clone();
    spanning.push(b'N');
    spanning.extend_from_slice(&right);

    let (across, _) = sketch_sequence(&spanning, params(20));
    let (whole, _) = sketch_sequence(&[left.clone(), right.clone()].concat(), params(20));
    let mut apart = [sketch_sequence(&left, params(20)).0, sketch_sequence(&right, params(20)).0].concat();
    apart.sort_unstable();
    apart.dedup();

    assert_eq!(across, apart);
    assert_ne!(across, whole);
}
