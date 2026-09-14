#![allow(dead_code)]

use std::collections::HashSet;

use rosella::quality::{Quality, Scorer};

/// Worth orders on bases, and every bin lands between the ladder's last rung and its first, so
/// the pool mechanics these tests cover run without a marker annotation behind them.
pub struct BasesScorer {
    lengths: Vec<usize>,
    span: usize,
}

const FLOOR: f64 = 55.0;
const BAND: f64 = 30.0;

impl BasesScorer {
    pub fn new(lengths: Vec<usize>) -> Self {
        let span = lengths.iter().sum();
        Self { lengths, span }
    }

    fn bases(&self, contigs: &[usize]) -> usize {
        contigs.iter().map(|contig| self.lengths[*contig]).sum()
    }
}

impl Scorer for BasesScorer {
    fn score(&self, contigs: &[usize]) -> Quality {
        Quality {
            completeness: FLOOR + BAND * self.bases(contigs) as f64 / self.span as f64,
            contamination: 0.0,
        }
    }

    fn features(&self, contigs: &[usize]) -> HashSet<u32> {
        contigs.iter().map(|contig| *contig as u32).collect()
    }

    fn sees_scale(&self) -> bool {
        false
    }
}
