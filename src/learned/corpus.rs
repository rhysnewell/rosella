use std::collections::HashMap;

use anyhow::Result;
use log::debug;
use ndarray::{Array2, ArrayView2, Axis, s};
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::learned::device::MOST_MATES;
use crate::learned::gemm::fresh;
use crate::learned::loss::Mates;
use crate::learned::neighbours::{CROSS, Lender};
use crate::learned::pieces::{Pieces, VIEWS, pieces, views};
use crate::learned::windows::SpanDepth;

const PSEUDO: f64 = 0.1;
const SPREAD_FLOOR: f64 = 1e-6;
const NO_BIN: u32 = u32::MAX;

pub struct Inputs<'a> {
    pub assembly: &'a str,
    pub names: &'a [String],
    pub lengths: &'a [usize],
    pub composition: ArrayView2<'a, f64>,
    pub coverage: ArrayView2<'a, f64>,
    pub kmer_size: usize,
    pub seed: u64,
    pub epochs: usize,
    pub patience: usize,
    pub stretches: usize,
    pub fast: bool,
    pub spans: &'a SpanDepth,
    pub families: Vec<Vec<u32>>,
    pub bins: Vec<Vec<u32>>,
    pub radius: bool,
    pub borrow: bool,
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

    pub(crate) fn blank(composition: usize, depth: usize) -> Self {
        let blank = |columns| Scale {
            centre: vec![0.0; columns],
            spread: vec![1.0; columns],
        };
        Self {
            composition: blank(composition),
            depth: blank(depth),
        }
    }

    pub(crate) fn parts(&self) -> [&[f32]; 4] {
        [
            &self.composition.centre,
            &self.composition.spread,
            &self.depth.centre,
            &self.depth.spread,
        ]
    }

    pub(crate) fn parts_mut(&mut self) -> [&mut [f32]; 4] {
        [
            &mut self.composition.centre,
            &mut self.composition.spread,
            &mut self.depth.centre,
            &mut self.depth.spread,
        ]
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

    fn depth_into(&self, depth: &[f32], out: &mut [f32]) {
        let logged = depth.iter().map(|depth| (f64::from(*depth) + PSEUDO).ln());
        self.depth.put(logged, out);
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
    bin_of: Vec<u32>,
    reach: Vec<f32>,
    lender: Option<Lender>,
}

pub struct Batch {
    pub features: Array2<f32>,
    pub labels: Vec<u32>,
    pub mates: Mates,
}

impl<'a> Corpus<'a> {
    pub fn new(inputs: &Inputs<'a>) -> Result<Self> {
        let contigs = inputs.composition.nrows();
        let wanted = (0..contigs)
            .map(|row| (inputs.names[row].as_str(), row))
            .collect::<HashMap<_, _>>();
        let crops = {
            let _timer = crate::timing::scope("learned_pieces");
            pieces(
                inputs.assembly,
                &wanted,
                inputs.kmer_size,
                inputs.seed,
                VIEWS,
                views,
            )?
        };
        let scales = Scales::fit(inputs.composition, inputs.coverage);
        let reach = match inputs.radius {
            true => reach(&scales, inputs.spans, contigs),
            false => Vec::new(),
        };
        let mut bin_of = vec![NO_BIN; if inputs.bins.is_empty() { 0 } else { contigs }];
        for (bin, members) in inputs.bins.iter().enumerate() {
            members
                .iter()
                .for_each(|row| bin_of[*row as usize] = bin as u32);
        }
        if !inputs.bins.is_empty() {
            let held = inputs.bins.iter().map(Vec::len).sum::<usize>();
            let bins = inputs.bins.len();
            debug!("Bin positives: {held} of {contigs} contigs in {bins} bins");
        }
        if !reach.is_empty() {
            let mut sorted = reach.clone();
            sorted.sort_unstable_by(f32::total_cmp);
            debug!(
                "Depth radius: median {:.2} z",
                sorted[sorted.len() / 2].sqrt()
            );
        }
        Ok(Self {
            scales,
            contigs,
            crops,
            spans: inputs.spans,
            width: inputs.composition.ncols(),
            bin_of,
            reach,
            lender: inputs.borrow.then(|| lender(inputs)),
        })
    }

    pub fn batch(&self, rng: &mut StdRng, contigs: &[u32], stretches: usize, near: bool) -> Batch {
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
        let key = self.lender.as_ref().map(|_| rng.random::<u64>());
        let columns = self.scales.width();
        let mut features = Array2::zeros((picks.len(), columns));
        features
            .as_slice_mut()
            .expect("a fresh table is contiguous")
            .par_chunks_mut(columns)
            .zip(&picks)
            .enumerate()
            .for_each(|(at, (out, (row, view)))| {
                let slot = row * VIEWS + view;
                let composition =
                    &self.crops.composition[slot * self.width..(slot + 1) * self.width];
                let mut depth = self
                    .spans
                    .view(*row, *view)
                    .iter()
                    .map(|value| f64::from(*value))
                    .collect::<Vec<_>>();
                if let (Some(lender), Some(key)) = (&self.lender, key) {
                    let mut rng = StdRng::seed_from_u64(crate::digest::mix(key ^ at as u64));
                    lender.lend(&mut rng, CROSS, *row, &mut depth);
                }
                self.scales.features(
                    composition.iter().map(|value| f64::from(*value)),
                    &depth,
                    out,
                );
            });
        let labels = picks.iter().map(|(row, _)| *row as u32).collect();
        let mates = self.mates(contigs, &features, 1 + stretches, near);
        Batch {
            features,
            labels,
            mates,
        }
    }

    // Bin mates past what the device loss holds are only left out of the negatives.
    fn mates(&self, contigs: &[u32], features: &Array2<f32>, views: usize, near: bool) -> Mates {
        let mut mates = Mates::default();
        let rows = |at: &usize| (at * views..(at + 1) * views).map(|row| row as u32);
        let bin = |at: usize| match self.bin_of.is_empty() {
            true => NO_BIN,
            false => self.bin_of[contigs[at] as usize],
        };
        let mut bins = HashMap::<u32, Vec<usize>>::new();
        for at in (0..contigs.len()).filter(|at| bin(*at) != NO_BIN) {
            bins.entry(bin(at)).or_default().push(at);
        }
        let most = MOST_MATES / views - 1;
        for members in bins.values().filter(|members| members.len() > 1) {
            for at in members {
                let others = members.iter().filter(|other| *other != at);
                let label = contigs[*at];
                mates
                    .pulled
                    .insert(label, others.clone().take(most).flat_map(rows).collect());
                mates
                    .masked
                    .insert(label, others.skip(most).flat_map(rows).collect());
            }
        }
        if !near || self.reach.is_empty() {
            return mates;
        }
        let whole = features
            .slice(s![..;views, self.scales.split()..])
            .to_owned();
        let turned = whole.t().as_standard_layout().into_owned();
        let dots = fresh(1.0, whole.view(), turned.view());
        let lengths = dots.diag().to_vec();
        let near = (0..contigs.len())
            .into_par_iter()
            .map(|at| {
                let reach = self.reach[contigs[at] as usize];
                let inside = (0..contigs.len())
                    .filter(|other| *other != at && (bin(at) == NO_BIN || bin(*other) != bin(at)))
                    .filter(|other| {
                        lengths[at] + lengths[*other] - 2.0 * dots[[at, *other]] <= reach
                    })
                    .flat_map(|other| rows(&other))
                    .collect::<Vec<_>>();
                (contigs[at], inside)
            })
            .collect::<Vec<_>>();
        for (label, inside) in near.into_iter().filter(|(_, inside)| !inside.is_empty()) {
            mates.masked.entry(label).or_default().extend(inside);
        }
        mates
    }
}

fn lender(inputs: &Inputs) -> Lender {
    let (contigs, samples) = (inputs.composition.nrows(), inputs.coverage.ncols() / 2);
    let unit = inputs.composition.mapv(|value| value as f32);
    let depth = (0..contigs)
        .flat_map(|row| depth_of(inputs.coverage, row))
        .collect();
    let labels = (0..contigs as u32).collect::<Vec<_>>();
    Lender::new(unit, &labels, depth, samples)
}

// A contig's own stretches are its depth noise, so any contig whose whole depth sits that close
// cannot be told from it by depth and is no negative (Huynh et al. 2022). Kept squared.
fn reach(scales: &Scales, spans: &SpanDepth, contigs: usize) -> Vec<f32> {
    let samples = scales.depth.centre.len();
    (0..contigs)
        .into_par_iter()
        .map(|row| {
            let (mut whole, mut view) = (vec![0.0; samples], vec![0.0; samples]);
            scales.depth_into(spans.view(row, 0), &mut whole);
            (1..VIEWS)
                .map(|at| {
                    scales.depth_into(spans.view(row, at), &mut view);
                    whole
                        .iter()
                        .zip(&view)
                        .map(|(whole, view)| (whole - view).powi(2))
                        .sum::<f32>()
                })
                .fold(0.0, f32::max)
        })
        .collect()
}

const EMBED_ROWS: usize = 1 << 16;

// Rows come back centred and unit length, so the composition distance 1 - rho reads 1 - cos.
pub fn embedded(
    scales: &Scales,
    width: usize,
    composition: ArrayView2<f64>,
    depth: impl Fn(usize) -> Vec<f64> + Sync,
    embed: impl Fn(ArrayView2<f32>) -> Array2<f32>,
) -> Array2<f64> {
    let columns = scales.width();
    let mut out = Array2::zeros((composition.nrows(), width));
    for (at, (chunk, mut target)) in composition
        .axis_chunks_iter(Axis(0), EMBED_ROWS)
        .zip(out.axis_chunks_iter_mut(Axis(0), EMBED_ROWS))
        .enumerate()
    {
        let mut input = Array2::zeros((chunk.nrows(), columns));
        input
            .as_slice_mut()
            .expect("a fresh table is contiguous")
            .par_chunks_mut(columns)
            .enumerate()
            .for_each(|(row, values)| {
                scales.features(
                    chunk.row(row).iter().copied(),
                    &depth(at * EMBED_ROWS + row),
                    values,
                )
            });
        let embedded = centred_unit(embed(input.view()));
        target.assign(&embedded.mapv(f64::from));
    }
    out
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
