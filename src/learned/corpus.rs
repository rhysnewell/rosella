use std::collections::HashMap;

use anyhow::Result;
use ndarray::{Array2, ArrayView2};
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::learned::neighbours::{NONE, RELATIVES, relatives};
use crate::learned::pieces::{Pieces, VIEWS, pieces};
use crate::learned::sample::poisson;
use crate::learned::windows::SpanDepth;

// A view's depth is read counts at its own span, so a short piece is as noisy as a short contig.
const READ: f64 = 150.0;
const PSEUDO: f64 = 0.1;
// Without window depth every view of a contig draws on one mean, so half the views take a share
// of a near genome's depth, as cross-mapped reads would lend it.
const MIXING: f64 = 0.5;
const SPREAD_FLOOR: f64 = 1e-6;

pub struct Inputs<'a> {
    pub assembly: &'a str,
    pub names: &'a [String],
    pub lengths: &'a [usize],
    pub composition: ArrayView2<'a, f64>,
    pub coverage: ArrayView2<'a, f64>,
    pub kmer_size: usize,
    pub seed: u64,
    pub spans: Option<&'a SpanDepth>,
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

fn resampled(rng: &mut impl Rng, depth: &[f64], span: f64) -> Vec<f64> {
    depth
        .iter()
        .map(|depth| poisson(rng, depth * span / READ) * READ / span)
        .collect()
}

pub struct Corpus<'a> {
    pub scales: Scales,
    pub contigs: usize,
    crops: Pieces,
    depth: Vec<f64>,
    peaks: Vec<f64>,
    relatives: Vec<[u32; RELATIVES]>,
    spans: Option<&'a SpanDepth>,
    width: usize,
    samples: usize,
}

pub struct Batch {
    pub features: Array2<f32>,
    pub labels: Vec<u32>,
}

impl<'a> Corpus<'a> {
    pub fn new(inputs: &Inputs<'a>) -> Result<Self> {
        let contigs = inputs.composition.nrows();
        let width = inputs.composition.ncols();
        let samples = inputs.coverage.ncols() / 2;
        let wanted = (0..contigs)
            .map(|row| (inputs.names[row].as_str(), row))
            .collect::<HashMap<_, _>>();
        let crops = {
            let _timer = crate::timing::scope("learned_pieces");
            pieces(inputs.assembly, &wanted, inputs.kmer_size, inputs.seed)?
        };
        let depth = (0..contigs)
            .flat_map(|row| depth_of(inputs.coverage, row))
            .collect::<Vec<_>>();
        let peaks = depth
            .chunks(samples.max(1))
            .map(|depth| depth.iter().copied().fold(1e-6, f64::max))
            .collect();
        let relatives = match inputs.spans {
            Some(_) => Vec::new(),
            None => {
                let _timer = crate::timing::scope("learned_relatives");
                let unit = inputs.composition.mapv(|value| value as f32);
                let labels = (0..contigs as u32).collect::<Vec<_>>();
                relatives(centred_unit(unit).view(), &labels)
            }
        };
        let scales = Scales {
            composition: Scale::fit(width, contigs, |row, column| {
                inputs.composition[[row, column]]
            }),
            depth: Scale::fit(samples, contigs, |row, sample| {
                (inputs.coverage[[row, 2 * sample]] + PSEUDO).ln()
            }),
        };
        Ok(Self {
            scales,
            contigs,
            crops,
            depth,
            peaks,
            relatives,
            spans: inputs.spans,
            width,
            samples,
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
        let key = rng.random::<u64>();
        let columns = self.scales.width();
        let mut features = Array2::zeros((picks.len(), columns));
        features
            .as_slice_mut()
            .expect("a fresh table is contiguous")
            .par_chunks_mut(columns)
            .zip(&picks)
            .enumerate()
            .for_each(|(at, (out, (row, view)))| {
                let mut rng = StdRng::seed_from_u64(crate::digest::mix(key ^ at as u64));
                let slot = row * VIEWS + view;
                let composition =
                    &self.crops.composition[slot * self.width..(slot + 1) * self.width];
                let span = f64::from(self.crops.spans[slot]);
                let depth = self.depth_view(&mut rng, *row, *view, span);
                self.scales.features(
                    composition.iter().map(|value| f64::from(*value)),
                    &depth,
                    out,
                );
            });
        let labels = picks.iter().map(|(row, _)| *row as u32).collect();
        Batch { features, labels }
    }

    fn depth_view(&self, rng: &mut StdRng, row: usize, view: usize, span: f64) -> Vec<f64> {
        if let Some(spans) = self.spans {
            return spans
                .view(row, view)
                .iter()
                .map(|value| f64::from(*value))
                .collect();
        }
        let samples = self.samples;
        let own = &self.depth[row * samples..(row + 1) * samples];
        let mut depth = own.to_vec();
        if rng.random::<f64>() < MIXING {
            let relative = self.relatives[row][rng.random_range(0..RELATIVES)];
            let share = rng.random::<f64>();
            if relative != NONE {
                let relative = relative as usize;
                let lent = &self.depth[relative * samples..(relative + 1) * samples];
                let scale = share * self.peaks[row] / self.peaks[relative];
                depth
                    .iter_mut()
                    .zip(lent)
                    .for_each(|(depth, lent)| *depth += scale * lent);
            }
        }
        resampled(rng, &depth, span)
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
