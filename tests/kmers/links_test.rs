use std::io::Write;

use rosella::kmers::links::{LinkScope, LinkSettings, links};
use rosella::kmers::sketch::{ContigSketches, SketchParams};

fn grow(seed: u64, length: usize) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 33) as usize % 4]
        })
        .collect()
}

fn drift(sequence: &[u8], every: usize) -> Vec<u8> {
    sequence
        .iter()
        .enumerate()
        .map(|(at, base)| {
            if at % every == 0 {
                b'A' + (base % 3)
            } else {
                *base
            }
        })
        .collect()
}

fn sketched(contigs: &[(&str, Vec<u8>)]) -> (tempfile::TempDir, ContigSketches) {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("assembly.fna");
    let mut file = std::fs::File::create(&path).unwrap();
    for (name, sequence) in contigs {
        writeln!(file, ">{name}").unwrap();
        file.write_all(sequence).unwrap();
        writeln!(file).unwrap();
    }
    drop(file);
    let sketches = ContigSketches::build(
        path.to_str().unwrap(),
        SketchParams {
            kmer_size: 31,
            scale: 1,
        },
    )
    .unwrap();
    (directory, sketches)
}

#[test]
fn near_identical_over_both_lengths_is_held_apart() {
    let one = grow(11, 4000);
    let other = drift(&one, 800);
    let (_directory, sketches) = sketched(&[("one", one), ("other", other)]);

    let found = links(&sketches, LinkSettings::default());

    assert_eq!(found.apart, vec![(0, 1)]);
    assert!(found.together.is_empty());
}

#[test]
fn a_fragment_of_one_longer_contig_is_held_together() {
    let host = grow(7, 4000);
    let inside = host[500..1200].to_vec();
    let (_directory, sketches) = sketched(&[("host", host), ("inside", inside)]);

    let found = links(&sketches, LinkSettings::default());

    assert!(found.apart.is_empty());
    assert_eq!(found.together, vec![(0, 1)]);
}

#[test]
fn a_fragment_two_contigs_both_hold_is_a_repeat_and_links_nothing() {
    let repeat = grow(3, 700);
    let mut one = grow(5, 4000);
    let mut other = grow(9, 4000);
    one.splice(500..500, repeat.iter().copied());
    other.splice(2000..2000, repeat.iter().copied());
    let (_directory, sketches) = sketched(&[("one", one), ("other", other), ("repeat", repeat)]);

    let found = links(&sketches, LinkSettings::default());

    assert!(found.apart.is_empty());
    assert!(found.together.is_empty());
}

#[test]
fn a_must_link_chain_collapses_to_one_component() {
    let host = grow(21, 6000);
    let head = host[200..1500].to_vec();
    let (_directory, sketches) =
        sketched(&[("host", host), ("head", head), ("apart", grow(31, 3000))]);

    let found = links(&sketches, LinkSettings::default());
    let components = found.components(3);

    assert_eq!(components[0], components[1]);
    assert_ne!(components[0], components[2]);
}

#[test]
fn a_scope_of_one_half_leaves_the_other_empty() {
    let one = grow(11, 4000);
    let other = drift(&one, 800);
    let contigs = [("one", one), ("other", other)];
    let (_directory, sketches) = sketched(&contigs);

    let held = links(
        &sketches,
        LinkSettings {
            scope: LinkScope::Together,
            ..LinkSettings::default()
        },
    );

    assert!(held.apart.is_empty());
}
