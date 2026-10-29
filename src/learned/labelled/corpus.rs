use std::collections::HashMap;

use anyhow::Result;
use ndarray::{Array2, ArrayView2};
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::kmers::halves::Halves;
use crate::learned::corpus::{Batch, Scales, centred_unit, depth_of};
use crate::learned::labelled::pieces::{CUTS, SIDES, pieces, sides};
use crate::learned::neighbours::{NONE, RELATIVES, relatives};
use crate::learned::pieces::{Pieces, view_key};
use crate::learned::sample::poisson;
use crate::learned::windows::{SpanDepth, Windows};

const OPTIONS: usize = 3;
// A view's depth is read counts at its own span, so a short piece is as noisy as a short contig.
const READ: f64 = 150.0;
const CROSS: f64 = 0.5;
const CROSS_SCALE: f64 = 1.0;

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
    pub epochs: usize,
    pub patience: usize,
    pub fast: bool,
    pub windows: Option<&'a Windows>,
}

pub fn resampled(rng: &mut impl Rng, depth: &[f64], span: f64) -> Vec<f64> {
    depth
        .iter()
        .map(|depth| poisson(rng, depth * span / READ) * READ / span)
        .collect()
}

// Each half reads its own depth at half the span, so a pair is as noisy as two contigs that long
// and the weight fit does not hand the space a free match.
pub fn half_depths(coverage: &[&[f64]], lengths: &[usize], seed: u64, side: u64) -> Vec<Vec<f64>> {
    coverage
        .iter()
        .zip(lengths)
        .enumerate()
        .map(|(at, (coverage, length))| {
            let key = crate::digest::mix(seed ^ (2 * at as u64 + side));
            let own = coverage.iter().step_by(2).copied().collect::<Vec<_>>();
            resampled(&mut StdRng::seed_from_u64(key), &own, *length as f64 / 2.0)
        })
        .collect()
}

pub struct Corpus {
    pub scales: Scales,
    pub held_rows: usize,
    pub others: Vec<usize>,
    pub held: Vec<usize>,
    labels: Vec<u32>,
    members: Vec<Vec<u32>>,
    options: Vec<f32>,
    option_spans: Vec<f32>,
    crops: Pieces,
    depth: Vec<f64>,
    peaks: Vec<f64>,
    spans: Option<SpanDepth>,
    relatives: Vec<[u32; RELATIVES]>,
    width: usize,
    samples: usize,
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

        let spans = match inputs.windows {
            Some(windows) => {
                let _timer = crate::timing::scope("learned_windows");
                let names = crate::rows::reorder(inputs.names, &rows);
                Some(windows.cut(&names, |at| {
                    let length = inputs.lengths[rows[at]];
                    match at < held.len() {
                        true => option_cuts(length, length >= inputs.halves_floor),
                        false => sides(view_key(inputs.seed, &names[at]), length),
                    }
                })?)
            }
            None => None,
        };
        let depth = rows
            .iter()
            .flat_map(|row| depth_of(inputs.coverage, *row))
            .collect::<Vec<_>>();
        let peaks = depth
            .chunks(samples.max(1))
            .map(|depth| depth.iter().copied().fold(1e-6, f64::max))
            .collect();
        let scales = Scales::fit(inputs.composition, inputs.coverage);
        let unit = Array2::from_shape_fn((rows.len(), width), |(at, column)| {
            inputs.composition[[rows[at], column]] as f32
        });
        let relatives = {
            let _timer = crate::timing::scope("learned_relatives");
            relatives(centred_unit(unit).view(), &labels)
        };
        Ok(Self {
            scales,
            held_rows: held.len(),
            others,
            held,
            labels,
            members,
            options,
            option_spans,
            crops,
            depth,
            peaks,
            spans,
            relatives,
            width,
            samples,
        })
    }

    pub fn labelled_bins(&self) -> usize {
        self.members.len()
    }

    pub fn bin_rows(&self, per_bin: usize) -> usize {
        self.members
            .iter()
            .map(|members| members.len().min(per_bin))
            .sum()
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
            let row = self.held_rows + *other as usize;
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
                let depth = self.depth_view(&mut rng, *row, *option, span);
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
                let slot = (row - self.held_rows) * SIDES + option;
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

    // Window depth already carries the read noise of the view's span, so it is not drawn again.
    fn depth_view(&self, rng: &mut StdRng, row: usize, option: usize, span: f64) -> Vec<f64> {
        let samples = self.samples;
        let mut depth = match &self.spans {
            Some(spans) => spans
                .view(row, option)
                .iter()
                .map(|value| f64::from(*value))
                .collect(),
            None => self.depth[row * samples..(row + 1) * samples].to_vec(),
        };
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
        match self.spans {
            Some(_) => depth,
            None => resampled(rng, &depth, span),
        }
    }
}

fn option_cuts(length: usize, halved: bool) -> Vec<(usize, usize)> {
    match halved {
        true => vec![(0, length), (0, length / 2), (length / 2, length)],
        false => vec![(0, length); OPTIONS],
    }
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
    let [first, second] = inputs.halves.composition(&names)?;
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
