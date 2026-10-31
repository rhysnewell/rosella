use std::collections::BTreeMap;
use std::io::Write;

use rosella::kmers::scan::{Floors, scan};
use rosella::refine::copies::copies;

const LONG: usize = 20_000;

fn random(seed: u64, length: usize) -> Vec<u8> {
    let mut state = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 62) as usize]
        })
        .collect()
}

fn rotated(sequence: &[u8], by: usize) -> Vec<u8> {
    sequence[by..]
        .iter()
        .chain(&sequence[..by])
        .copied()
        .collect()
}

#[test]
fn only_contigs_repeated_inside_their_own_bin_leave() {
    let element = random(1, LONG);
    let sequences = [
        element.clone(),
        rotated(&element, 7_000),
        random(2, LONG),
        element[..5_000].to_vec(),
        rotated(&element, 13_000),
        random(3, LONG),
    ];
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("assembly.fasta");
    let mut sink = std::fs::File::create(&path).unwrap();
    let names = (0..sequences.len())
        .map(|at| format!("contig_{at}"))
        .collect::<Vec<_>>();
    for (name, sequence) in names.iter().zip(&sequences) {
        writeln!(sink, ">{name}").unwrap();
        sink.write_all(sequence).unwrap();
        writeln!(sink).unwrap();
    }
    drop(sink);
    let floors = Floors {
        sketch: Some(0),
        ..Floors::default()
    };
    let mut sketches = scan(&path.to_string_lossy(), 4, floors)
        .unwrap()
        .sketches
        .unwrap();
    sketches.align_to(&names).unwrap();
    let lengths = sequences.iter().map(Vec::len).collect::<Vec<_>>();
    let bins = BTreeMap::from([(0, vec![0, 1, 2, 3]), (1, vec![4, 5])]);

    assert_eq!(copies(&bins, &sketches, &lengths, 10_000), vec![0, 1]);
}
