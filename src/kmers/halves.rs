use std::collections::HashMap;

use anyhow::{Result, anyhow};
use ndarray::Array2;

use crate::kmers::kmer_counting::frequencies;

// Counts take half the memory of frequencies. They give the same rows to the bit.
#[derive(Default)]
pub struct Halves {
    width: usize,
    index: HashMap<String, usize>,
    lengths: Vec<usize>,
    counts: Vec<u32>,
    kmers: Vec<u32>,
}

pub(crate) type Counted = (usize, [(Vec<u32>, u32); 2]);

impl Halves {
    pub(crate) fn new(width: usize) -> Self {
        Self {
            width,
            ..Self::default()
        }
    }

    pub(crate) fn push(&mut self, name: String, (length, halves): Counted) {
        self.index.insert(name, self.lengths.len());
        self.lengths.push(length);
        for (counts, kmers) in halves {
            self.counts.extend(counts);
            self.kmers.push(kmers);
        }
    }

    pub fn composition(&self, names: &[&str], kmer_size: usize) -> Result<[Array2<f64>; 2]> {
        let rows = self.rows(names)?;
        Ok([
            self.table(&rows, &[0], kmer_size)?,
            self.table(&rows, &[1], kmer_size)?,
        ])
    }

    // The k-mer table holds learned rows once the space trains, so a test of k-mers alone reads
    // them from here, short only the k-mers that span the cut.
    pub fn whole(&self, names: &[&str], kmer_size: usize) -> Result<Array2<f64>> {
        self.table(&self.rows(names)?, &[0, 1], kmer_size)
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

    fn table(&self, rows: &[usize], halves: &[usize], kmer_size: usize) -> Result<Array2<f64>> {
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
        crate::kmers::clr::clr(&mut table, &lengths, kmer_size)?;
        Ok(table)
    }
}
