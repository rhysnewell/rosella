#![allow(dead_code)]
//! A contig set with real sequence behind it, because `ContigSketches::build` reads a file.

use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::kmers::sketch::{ContigSketches, SketchParams};

pub const SCALE: u64 = 20;
static FIXTURE: AtomicUsize = AtomicUsize::new(0);

pub fn grow(seed: u64, length: usize) -> Vec<u8> {
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

/// A sibling strain: the same sequence with a substitution every hundredth base.
pub fn sibling(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .enumerate()
        .map(|(at, base)| {
            if at % 100 == 0 {
                b'A' + (base % 3)
            } else {
                *base
            }
        })
        .collect()
}

pub struct Fixture {
    pub sketches: ContigSketches,
    pub lengths: Vec<usize>,
    pub coverage: Array2<f64>,
    pub tnf: Array2<f64>,
}

impl Fixture {
    pub fn new(name: &str, contigs: Vec<Vec<u8>>) -> Self {
        let unique = FIXTURE.fetch_add(1, Ordering::Relaxed);
        let path =
            std::env::temp_dir().join(format!("rosella_{name}_{}_{unique}.fa", std::process::id()));
        let mut handle = std::fs::File::create(&path).expect("temp fasta");
        for (index, sequence) in contigs.iter().enumerate() {
            writeln!(handle, ">c{index}").unwrap();
            handle.write_all(sequence).unwrap();
            writeln!(handle).unwrap();
        }
        drop(handle);
        let sketches = ContigSketches::build(
            path.to_str().unwrap(),
            SketchParams {
                kmer_size: 31,
                scale: SCALE,
            },
        )
        .expect("sketches");
        std::fs::remove_file(&path).ok();
        let lengths = contigs.iter().map(Vec::len).collect::<Vec<_>>();
        Self {
            coverage: Array2::zeros((lengths.len(), 2)),
            tnf: Array2::zeros((lengths.len(), 2)),
            sketches,
            lengths,
        }
    }

    pub fn features(&self) -> ContigFeatures<'_> {
        ContigFeatures::new(&self.coverage, &self.tnf, &self.lengths)
            .with_sketches(Some(&self.sketches))
    }
}
