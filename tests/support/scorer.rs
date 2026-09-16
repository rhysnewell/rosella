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
            ..Default::default()
        }
    }

    fn features(&self, contigs: &[usize]) -> HashSet<u32> {
        contigs.iter().map(|contig| *contig as u32).collect()
    }
}

/// Completeness and contamination against a known contig to genome map, so a split of a fused
/// bin can outscore the bin, which is the case worth ordering has to get right.
pub struct GenomeScorer {
    genomes: Vec<usize>,
}

impl GenomeScorer {
    pub fn new(genomes: Vec<usize>) -> Self {
        Self { genomes }
    }

    fn held(&self, genome: usize) -> usize {
        self.genomes.iter().filter(|held| **held == genome).count()
    }
}

impl Scorer for GenomeScorer {
    fn score(&self, contigs: &[usize]) -> Quality {
        let mut counts = std::collections::BTreeMap::new();
        for contig in contigs {
            *counts.entry(self.genomes[*contig]).or_insert(0usize) += 1;
        }
        let Some((genome, held)) = counts.into_iter().max_by_key(|(_, held)| *held) else {
            return Quality::default();
        };
        Quality {
            completeness: 100.0 * held as f64 / self.held(genome) as f64,
            contamination: 100.0 * (contigs.len() - held) as f64 / held as f64,
            ..Default::default()
        }
    }

    fn features(&self, contigs: &[usize]) -> HashSet<u32> {
        contigs.iter().map(|contig| *contig as u32).collect()
    }
}

/// Marker families by contig, so a test can make two pieces share a family or not. The peel
/// gate turns on that overlap and nothing else does.
pub struct FamilyScorer {
    families: Vec<Option<u32>>,
}

impl FamilyScorer {
    pub fn new(families: Vec<Option<u32>>) -> Self {
        Self { families }
    }
}

impl Scorer for FamilyScorer {
    fn score(&self, contigs: &[usize]) -> Quality {
        Quality {
            completeness: self.features(contigs).len() as f64,
            ..Default::default()
        }
    }

    fn features(&self, contigs: &[usize]) -> HashSet<u32> {
        contigs
            .iter()
            .filter_map(|contig| self.families[*contig])
            .collect()
    }
}
