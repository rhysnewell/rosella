//! Gene calling, which the gene family scorer reads before it can score anything.

use rosella::quality::orfs::{translate, translate_reverse};

#[test]
fn translation_rewrites_the_start_and_drops_the_stop() {
    assert_eq!(translate(b"GTGAAATAA", true), "MK");
    assert_eq!(translate(b"GTGAAATAA", false), "VK");
    assert_eq!(translate(b"ATGNNNTAA", true), "MX");
}

#[test]
fn the_reverse_strand_matches_translating_the_reverse_complement() {
    fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|byte| match byte.to_ascii_uppercase() {
                b'A' => b'T',
                b'T' | b'U' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                other => other,
            })
            .collect()
    }

    let alphabet = b"ACGTNacgtuU";
    let mut state = 0x2545f4914f6cdd1du64;
    for length in 1..96usize {
        let coding = (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                alphabet[(state % alphabet.len() as u64) as usize]
            })
            .collect::<Vec<u8>>();
        for complete_start in [true, false] {
            assert_eq!(
                translate_reverse(&coding, complete_start),
                translate(&reverse_complement(&coding), complete_start),
                "length {length} start {complete_start}"
            );
        }
    }
}
