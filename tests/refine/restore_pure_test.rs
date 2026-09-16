//! The k-mer band that keeps a bin the marker table condemns for its own duplicated families.

use std::collections::HashSet;
use std::io::Write;

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::kmers::sketch::ContigSketches;
use rosella::quality::{Quality, Scorer};
use rosella::refine::restore::{Judge, restore};
use rosella::refine::rung::Rung;

const PIECE: usize = 200_000;
const REPEAT: usize = 36_000;
const CONTIGS: usize = 10;

/// An evolved strain reads complete and contaminated whole and in pieces alike, so worth alone
/// can never separate the two arrangements and only the k-mers can.
struct EvoScorer;

impl Scorer for EvoScorer {
    fn score(&self, _contigs: &[usize]) -> Quality {
        Quality {
            completeness: 100.0,
            contamination: 8.33,
        }
    }

    fn features(&self, contigs: &[usize]) -> HashSet<u32> {
        contigs.iter().map(|contig| *contig as u32).collect()
    }
}

fn accept() -> Rung {
    Rung {
        floor: 1,
        completeness: 90.0,
        contamination: 5.0,
    }
}

fn reported() -> Rung {
    Rung {
        floor: 1,
        completeness: 50.0,
        contamination: f64::INFINITY,
    }
}

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
        .map(|(at, base)| match at % every {
            0 => b'A' + (base % 3),
            _ => *base,
        })
        .collect()
}

/// Contigs 0 to 3 are one genome repeating 4.5 per cent of itself, the rate a duplicated marker
/// family gives; 4 to 7 are unrelated genomes sharing nothing; 8 and 9 are two close strains.
fn assembly(directory: &std::path::Path) -> (String, Vec<String>) {
    let path = directory.join("assembly.fna");
    let mut sink = std::fs::File::create(&path).unwrap();
    let first = grow(11, PIECE);
    let mut pieces = vec![first.clone(), grow(12, PIECE), grow(13, PIECE)];
    let mut repeated = grow(14, PIECE - REPEAT);
    repeated.extend_from_slice(&first[..REPEAT]);
    pieces.push(repeated);
    pieces.extend((20..24).map(|seed| grow(seed, PIECE)));
    let twin = grow(31, PIECE);
    pieces.push(drift(&twin, 100));
    pieces.push(twin);

    let mut names = Vec::new();
    for (at, sequence) in pieces.iter().enumerate() {
        let name = format!("contig_{at}");
        writeln!(sink, ">{name}").unwrap();
        sink.write_all(sequence).unwrap();
        writeln!(sink).unwrap();
        names.push(name);
    }
    (path.to_string_lossy().into_owned(), names)
}

macro_rules! held {
    ($dissolved:expr, $promoted:expr) => {{
        let directory = tempfile::tempdir().unwrap();
        let (path, names) = assembly(directory.path());
        let mut sketches = ContigSketches::build(&path).unwrap();
        sketches.align_to(&names).unwrap();

        let coverage = Array2::zeros((CONTIGS, 2));
        let tnf = Array2::zeros((CONTIGS, 2));
        let lengths = vec![PIECE; CONTIGS];
        let features =
            ContigFeatures::new(&coverage, &tnf, &lengths).with_sketches(Some(&sketches));
        let held = Judge {
            features: &features,
            quality: &EvoScorer,
            reported: reported(),
            accept: accept(),
        };
        restore(&held, 2.0, &$dissolved, $promoted)
    }};
}

#[test]
fn one_genome_repeating_its_own_families_goes_back_whole() {
    let held = held!(
        vec![(0usize, vec![0, 1, 2, 3])],
        vec![vec![0, 1], vec![2, 3]]
    );

    assert_eq!(held.bins, 1);
    assert!(held.promoted.is_empty());
}

/// Unrelated genomes in one bin repeat nothing either, so duplication alone would keep them.
/// The floor is what tells them apart from a strain that duplicates its own families.
#[test]
fn a_bin_of_unrelated_genomes_falls_under_the_floor() {
    let promoted = vec![vec![4, 5], vec![6, 7]];
    let held = held!(vec![(0usize, vec![4, 5, 6, 7])], promoted.clone());

    assert_eq!(held.bins, 0);
    assert_eq!(held.promoted, promoted);
}

/// Two close strains share far more than a genome repeats of itself, which puts them over the bar.
#[test]
fn a_bin_holding_two_strains_sits_over_the_bar() {
    let promoted = vec![vec![8], vec![9]];
    let held = held!(vec![(0usize, vec![8, 9])], promoted.clone());

    assert_eq!(held.bins, 0);
    assert_eq!(held.promoted, promoted);
}
