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
        let rows = names
            .iter()
            .map(|name| {
                self.index
                    .get(*name)
                    .copied()
                    .ok_or_else(|| anyhow!("{name} was not counted in halves"))
            })
            .collect::<Result<Vec<_>>>()?;
        let side = |half: usize| -> Result<Array2<f64>> {
            let mut values = Vec::with_capacity(rows.len() * self.width);
            let mut lengths = Vec::with_capacity(rows.len());
            for row in &rows {
                let slot = 2 * row + half;
                let counts = &self.counts[slot * self.width..(slot + 1) * self.width];
                values.extend(frequencies(counts, self.kmers[slot]));
                let length = self.lengths[*row];
                lengths.push(if half == 0 {
                    length / 2
                } else {
                    length - length / 2
                });
            }
            let mut table = Array2::from_shape_vec((rows.len(), self.width), values)?;
            crate::kmers::clr::clr(&mut table, &lengths, kmer_size)?;
            Ok(table)
        };
        Ok([side(0)?, side(1)?])
    }
}
