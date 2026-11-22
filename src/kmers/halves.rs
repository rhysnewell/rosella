use std::collections::HashMap;

use anyhow::{Result, anyhow};
use ndarray::Array2;

use crate::kmers::kmer_counting::{Block, frequencies};

// Counts take half the memory of frequencies. They give the same rows to the bit.
#[derive(Default)]
pub struct Halves {
    width: usize,
    kmer_size: usize,
    gc_bases: Vec<u32>,
    index: HashMap<String, usize>,
    lengths: Vec<usize>,
    counts: Vec<u32>,
    kmers: Vec<u32>,
}

pub(crate) type Counted = (usize, [(Vec<u32>, u32); 2]);

impl Halves {
    pub(crate) fn new(block: &Block) -> Self {
        Self {
            width: block.width,
            kmer_size: block.kmer_size(),
            gc_bases: block.gc_bases(),
            ..Self::default()
        }
    }

    // The bases are gone once the scan ends, and every base sits in as many k-mers as the
    // k-mer is long, so the counts give the GC share back.
    pub fn gc(&self, name: &str) -> Option<f64> {
        let row = *self.index.get(name)?;
        let side = &self.counts[2 * row * self.width..(2 * row + 2) * self.width];
        let strong = side
            .iter()
            .zip(self.gc_bases.iter().cycle())
            .map(|(count, bases)| u64::from(*count) * u64::from(*bases))
            .sum::<u64>();
        let kmers = u64::from(self.kmers[2 * row]) + u64::from(self.kmers[2 * row + 1]);
        (kmers > 0).then(|| strong as f64 / (kmers * self.kmer_size as u64) as f64)
    }

    pub(crate) fn push(&mut self, name: String, (length, halves): Counted) {
        self.index.insert(name, self.lengths.len());
        self.lengths.push(length);
        for (counts, kmers) in halves {
            self.counts.extend(counts);
            self.kmers.push(kmers);
        }
    }

    pub fn composition(&self, names: &[&str]) -> Result<[Array2<f64>; 2]> {
        let rows = self.rows(names)?;
        Ok([self.table(&rows, &[0])?, self.table(&rows, &[1])?])
    }

    // The k-mer table holds learned rows once the space trains, so a test of k-mers alone reads
    // them from here, short only the k-mers that span the cut.
    pub fn whole(&self, names: &[&str]) -> Result<Array2<f64>> {
        self.table(&self.rows(names)?, &[0, 1])
    }

    fn rows(&self, names: &[&str]) -> Result<Vec<usize>> {
        names
            .iter()
            .map(|name| {
                self.index
                    .get(*name)
                    .copied()
                    .ok_or_else(|| anyhow!("{name} was not counted in halves"))
            })
            .collect()
    }

    fn table(&self, rows: &[usize], halves: &[usize]) -> Result<Array2<f64>> {
        let mut values = Vec::with_capacity(rows.len() * self.width);
        let mut lengths = Vec::with_capacity(rows.len());
        let mut counts = vec![0; self.width];
        for row in rows {
            counts.fill(0);
            let mut kmers = 0;
            let mut length = 0;
            for half in halves {
                let slot = 2 * row + half;
                let side = &self.counts[slot * self.width..(slot + 1) * self.width];
                for (sum, count) in counts.iter_mut().zip(side) {
                    *sum += count;
                }
                kmers += self.kmers[slot];
                let whole = self.lengths[*row];
                length += if *half == 0 {
                    whole / 2
                } else {
                    whole - whole / 2
                };
            }
            values.extend(frequencies(&counts, kmers));
            lengths.push(length);
        }
        let mut table = Array2::from_shape_vec((rows.len(), self.width), values)?;
        crate::kmers::clr::clr(&mut table, &lengths, self.kmer_size)?;
        Ok(table)
    }
}
