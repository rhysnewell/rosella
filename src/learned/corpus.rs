use std::collections::HashMap;

use anyhow::Result;
use ndarray::{Array2, ArrayView2};
use rand::{Rng, rngs::StdRng};
use rayon::prelude::*;

use crate::learned::pieces::{Pieces, VIEWS, pieces};
use crate::learned::windows::SpanDepth;

const PSEUDO: f64 = 0.1;
const SPREAD_FLOOR: f64 = 1e-6;

pub struct Inputs<'a> {
    pub assembly: &'a str,
    pub names: &'a [String],
    pub lengths: &'a [usize],
    pub composition: ArrayView2<'a, f64>,
    pub coverage: ArrayView2<'a, f64>,
    pub kmer_size: usize,
    pub seed: u64,
    pub epochs: usize,
    pub stretches: usize,
    pub spans: &'a SpanDepth,
}

pub struct Scale {
    centre: Vec<f32>,
    spread: Vec<f32>,
}

impl Scale {
    fn fit(columns: usize, rows: usize, value: impl Fn(usize, usize) -> f64 + Sync) -> Self {
        let (centre, spread) = (0..columns)
            .into_par_iter()
            .map(|column| {
                let mean = (0..rows).map(|row| value(row, column)).sum::<f64>() / rows as f64;
                let variance = (0..rows)
                    .map(|row| (value(row, column) - mean).powi(2))
                    .sum::<f64>()
                    / rows as f64;
                (mean as f32, (variance.sqrt() + SPREAD_FLOOR) as f32)
            })
            .unzip();
        Self { centre, spread }
    }

    fn put(&self, values: impl Iterator<Item = f64>, out: &mut [f32]) {
        for (((out, value), centre), spread) in out
            .iter_mut()
            .zip(values)
            .zip(&self.centre)
            .zip(&self.spread)
        {
            *out = (value as f32 - centre) / spread;
        }
    }
}

pub struct Scales {
    composition: Scale,
    depth: Scale,
}

impl Scales {
    pub fn fit(composition: ArrayView2<f64>, coverage: ArrayView2<f64>) -> Self {
        let (contigs, width) = composition.dim();
        Self {
            composition: Scale::fit(width, contigs, |row, column| composition[[row, column]]),
            depth: Scale::fit(coverage.ncols() / 2, contigs, |row, sample| {
                (coverage[[row, 2 * sample]] + PSEUDO).ln()
            }),
        }
    }

    pub fn width(&self) -> usize {
        self.composition.centre.len() + self.depth.centre.len()
    }

    pub fn split(&self) -> usize {
        self.composition.centre.len()
    }

    pub fn features(&self, composition: impl Iterator<Item = f64>, depth: &[f64], out: &mut [f32]) {
        let (left, right) = out.split_at_mut(self.composition.centre.len());
        self.composition.put(composition, left);
        self.depth
            .put(depth.iter().map(|depth| (depth + PSEUDO).ln()), right);
    }
}

pub fn depth_of(coverage: ArrayView2<f64>, row: usize) -> Vec<f64> {
    (0..coverage.ncols() / 2)
        .map(|sample| coverage[[row, 2 * sample]])
        .collect()
}

pub struct Corpus<'a> {
    pub scales: Scales,
    pub contigs: usize,
    crops: Pieces,
    spans: &'a SpanDepth,
    width: usize,
}

pub struct Batch {
    pub features: Array2<f32>,
    pub labels: Vec<u32>,
}

impl<'a> Corpus<'a> {
    pub fn new(inputs: &Inputs<'a>) -> Result<Self> {
        let contigs = inputs.composition.nrows();
        let wanted = (0..contigs)
            .map(|row| (inputs.names[row].as_str(), row))
            .collect::<HashMap<_, _>>();
        let crops = {
            let _timer = crate::timing::scope("learned_pieces");
            pieces(inputs.assembly, &wanted, inputs.kmer_size, inputs.seed)?
        };
        Ok(Self {
            scales: Scales::fit(inputs.composition, inputs.coverage),
            contigs,
            crops,
            spans: inputs.spans,
            width: inputs.composition.ncols(),
        })
    }

    pub fn batch(&self, rng: &mut StdRng, contigs: &[u32], stretches: usize) -> Batch {
        let mut picks = Vec::with_capacity(contigs.len() * (1 + stretches));
        for contig in contigs {
            let row = *contig as usize;
            picks.push((row, 0));
            let mut slots = (1..VIEWS).collect::<Vec<_>>();
            for at in 0..stretches {
                let pick = rng.random_range(at..slots.len());
                slots.swap(at, pick);
                picks.push((row, slots[at]));
            }
        }
        // Every view reads its window depth, so this draw goes unread. It stays so a seed trains
        // the same space.
        rng.random::<u64>();
        let columns = self.scales.width();
        let mut features = Array2::zeros((picks.len(), columns));
        features
            .as_slice_mut()
            .expect("a fresh table is contiguous")
            .par_chunks_mut(columns)
            .zip(&picks)
            .for_each(|(out, (row, view))| {
                let slot = row * VIEWS + view;
                let composition =
                    &self.crops.composition[slot * self.width..(slot + 1) * self.width];
                let depth = self
                    .spans
                    .view(*row, *view)
                    .iter()
                    .map(|value| f64::from(*value))
                    .collect::<Vec<_>>();
                self.scales.features(
                    composition.iter().map(|value| f64::from(*value)),
                    &depth,
                    out,
                );
            });
        let labels = picks.iter().map(|(row, _)| *row as u32).collect();
        Batch { features, labels }
    }
}

pub fn centred_unit(mut table: Array2<f32>) -> Array2<f32> {
    let columns = table.ncols();
    table
        .as_slice_mut()
        .expect("a fresh table is contiguous")
        .par_chunks_mut(columns)
        .for_each(|values| {
            let mean = values.iter().sum::<f32>() / columns as f32;
            values.iter_mut().for_each(|value| *value -= mean);
            let length = values.iter().map(|value| value * value).sum::<f32>().sqrt();
            values
                .iter_mut()
                .for_each(|value| *value /= length.max(1e-12));
        });
    table
}
