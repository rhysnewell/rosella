use std::collections::HashMap;

use anyhow::Result;
use ndarray::{Array2, ArrayView2};
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::kmers::halves::Halves;
use crate::learned::neighbours::{NONE, RELATIVES, relatives};
use crate::learned::pieces::{CUTS, Pieces, SIDES, pieces};
use crate::learned::sample::poisson;

const OPTIONS: usize = 3;
// A view's depth is read counts at its own span, so a short piece is as noisy as a short contig.
const READ: f64 = 150.0;
const PSEUDO: f64 = 0.1;
const CROSS: f64 = 0.5;
const CROSS_SCALE: f64 = 1.0;
const SPREAD_FLOOR: f64 = 1e-6;

pub struct Inputs<'a> {
    pub assembly: &'a str,
    pub names: &'a [String],
    pub lengths: &'a [usize],
    pub composition: ArrayView2<'a, f64>,
    pub coverage: ArrayView2<'a, f64>,
    pub labels: &'a [Vec<usize>],
    pub halves: &'a Halves,
    pub halves_floor: usize,
    pub kmer_size: usize,
    pub seed: u64,
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
    pub composition: Scale,
    pub depth: Scale,
}

impl Scales {
    pub fn width(&self) -> usize {
        self.composition.centre.len() + self.depth.centre.len()
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

pub fn resampled(rng: &mut impl Rng, depth: &[f64], span: f64) -> Vec<f64> {
    depth
        .iter()
        .map(|depth| poisson(rng, depth * span / READ) * READ / span)
        .collect()
}

pub struct Corpus {
    pub scales: Scales,
    pub bins: usize,
    pub others: usize,
    pub held: Vec<usize>,
    labels: Vec<u32>,
    members: Vec<Vec<u32>>,
    options: Vec<f32>,
    option_spans: Vec<f32>,
    crops: Pieces,
    depth: Vec<f64>,
    peaks: Vec<f64>,
    relatives: Vec<[u32; RELATIVES]>,
    width: usize,
    samples: usize,
}

pub struct Batch {
    pub features: Array2<f32>,
    pub labels: Vec<u32>,
}

impl Corpus {
    pub fn new(inputs: &Inputs) -> Result<Self> {
        let long = inputs.composition.nrows();
        let width = inputs.composition.ncols();
        let samples = inputs.coverage.ncols() / 2;
        let mut bin_of = vec![None; long];
        let kept = inputs
            .labels
            .iter()
            .map(|members| {
                members
                    .iter()
                    .copied()
                    .filter(|row| *row < long)
                    .collect::<Vec<_>>()
            })
            .filter(|members| members.len() >= 2)
            .collect::<Vec<_>>();
        for (bin, members) in kept.iter().enumerate() {
            for row in members {
                bin_of[*row] = Some(bin);
            }
        }
        let held = kept.iter().flatten().copied().collect::<Vec<_>>();
        let others = (0..long)
            .filter(|row| bin_of[*row].is_none())
            .collect::<Vec<_>>();
        let rows = held.iter().chain(&others).copied().collect::<Vec<_>>();
        let labels = held
            .iter()
            .map(|row| bin_of[*row].expect("a held row has a bin") as u32)
            .chain((0..others.len()).map(|at| (kept.len() + at) as u32))
            .collect::<Vec<_>>();
        let mut members = vec![Vec::new(); kept.len()];
        for (at, label) in labels[..held.len()].iter().enumerate() {
            members[*label as usize].push(at as u32);
        }

        let (options, option_spans) = bin_options(inputs, &held)?;
        let wanted = others
            .iter()
            .enumerate()
            .map(|(slot, row)| (inputs.names[*row].as_str(), slot))
            .collect::<HashMap<_, _>>();
        let crops = {
            let _timer = crate::timing::scope("learned_pieces");
            pieces(inputs.assembly, &wanted, inputs.kmer_size, inputs.seed)?
        };

        let depth = rows
            .iter()
            .flat_map(|row| depth_of(inputs.coverage, *row))
            .collect::<Vec<_>>();
        let peaks = depth
            .chunks(samples.max(1))
            .map(|depth| depth.iter().copied().fold(1e-6, f64::max))
            .collect();
        let scales = Scales {
            composition: Scale::fit(width, long, |row, column| inputs.composition[[row, column]]),
            depth: Scale::fit(samples, long, |row, sample| {
                (inputs.coverage[[row, 2 * sample]] + PSEUDO).ln()
            }),
        };
        let unit = Array2::from_shape_fn((rows.len(), width), |(at, column)| {
            inputs.composition[[rows[at], column]] as f32
        });
        let relatives = {
            let _timer = crate::timing::scope("learned_relatives");
            relatives(centred_unit(unit).view(), &labels)
        };
        Ok(Self {
            scales,
            bins: held.len(),
            others: others.len(),
            held,
            labels,
            members,
            options,
            option_spans,
            crops,
            depth,
            peaks,
            relatives,
            width,
            samples,
        })
    }

    pub fn labelled_bins(&self) -> usize {
        self.members.len()
    }

    pub fn batch(&self, rng: &mut StdRng, others: &[u32], per_bin: usize) -> Batch {
        let mut picks = Vec::new();
        for members in &self.members {
            let mut pool = members.clone();
            let take = per_bin.min(pool.len());
            for at in 0..take {
                let swap = rng.random_range(at..pool.len());
                pool.swap(at, swap);
                let option = rng.random_range(0..OPTIONS);
                picks.push((pool[at] as usize, option, false));
            }
        }
        for other in others {
            let cut = rng.random_range(0..CUTS);
            let side = rng.random_range(0..2);
            let row = self.bins + *other as usize;
            picks.push((row, 2 * cut + side, true));
            picks.push((row, 2 * cut + 1 - side, true));
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
            .for_each(|(at, (out, (row, option, cropped)))| {
                let mut rng = StdRng::seed_from_u64(crate::digest::mix(key ^ at as u64));
                let (composition, span) = self.view(*row, *option, *cropped);
                let depth = self.depth_view(&mut rng, *row, span);
                self.scales.features(
                    composition.iter().map(|value| f64::from(*value)),
                    &depth,
                    out,
                );
            });
        let labels = picks.iter().map(|(row, _, _)| self.labels[*row]).collect();
        Batch { features, labels }
    }

    fn view(&self, row: usize, option: usize, cropped: bool) -> (&[f32], f64) {
        let width = self.width;
        match cropped {
            true => {
                let slot = (row - self.bins) * SIDES + option;
                (
                    &self.crops.composition[slot * width..(slot + 1) * width],
                    f64::from(self.crops.spans[slot]),
                )
            }
            false => {
                let slot = row * OPTIONS + option;
                (
                    &self.options[slot * width..(slot + 1) * width],
                    f64::from(self.option_spans[slot]),
                )
            }
        }
    }

    fn depth_view(&self, rng: &mut StdRng, row: usize, span: f64) -> Vec<f64> {
        let samples = self.samples;
        let mut depth = self.depth[row * samples..(row + 1) * samples].to_vec();
        if rng.random::<f64>() < CROSS {
            let relative = self.relatives[row][rng.random_range(0..RELATIVES)];
            let share = rng.random::<f64>() * CROSS_SCALE;
            if relative != NONE {
                let relative = relative as usize;
                let lent = &self.depth[relative * samples..(relative + 1) * samples];
                let scale = share * self.peaks[row] / self.peaks[relative];
                for (depth, lent) in depth.iter_mut().zip(lent) {
                    *depth += scale * lent;
                }
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

// A label contig offers its whole composition or either half, so bin mates meet at more than
// one length.
fn bin_options(inputs: &Inputs, held: &[usize]) -> Result<(Vec<f32>, Vec<f32>)> {
    let width = inputs.composition.ncols();
    let halved = held
        .iter()
        .copied()
        .filter(|row| inputs.lengths[*row] >= inputs.halves_floor)
        .collect::<Vec<_>>();
    let names = halved
        .iter()
        .map(|row| inputs.names[*row].as_str())
        .collect::<Vec<_>>();
    let [first, second] = inputs.halves.composition(&names, inputs.kmer_size)?;
    let half_at = halved
        .iter()
        .enumerate()
        .map(|(at, row)| (*row, at))
        .collect::<HashMap<_, _>>();
    let mut options = Vec::with_capacity(held.len() * OPTIONS * width);
    let mut spans = Vec::with_capacity(held.len() * OPTIONS);
    for row in held {
        let whole = inputs.composition.row(*row);
        let length = inputs.lengths[*row] as f32;
        match half_at.get(row) {
            Some(at) => {
                options.extend(whole.iter().map(|value| *value as f32));
                options.extend(first.row(*at).iter().map(|value| *value as f32));
                options.extend(second.row(*at).iter().map(|value| *value as f32));
                spans.extend([length, length / 2.0, length / 2.0]);
            }
            None => {
                for _ in 0..OPTIONS {
                    options.extend(whole.iter().map(|value| *value as f32));
                }
                spans.extend([length; OPTIONS]);
            }
        }
    }
    Ok((options, spans))
}
