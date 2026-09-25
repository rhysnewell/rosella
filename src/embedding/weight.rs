use rand::{SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::embedding::metrics::{MIN_VAR, combine, metabat_with, rho};

pub mod noise;

use noise::Partners;

pub const NEIGHBOURS: usize = 10;
const STEPS: usize = 20;
const STEP: f64 = 0.05;
const PLATEAU: f64 = 0.005;
const SLOPES: usize = 41;
const SLOPE_REACH: f64 = 0.5;

pub struct Contigs<'a> {
    pub coverage: Vec<&'a [f64]>,
    pub whole: Vec<&'a [f64]>,
    pub first: Vec<Vec<f64>>,
    pub second: Vec<Vec<f64>>,
    pub lengths: Vec<usize>,
}

pub struct Derived {
    pub weight: f64,
    pub line: Line,
}

// The weight at the reference length, moving by the slope per tenfold of the shorter contig.
#[derive(Debug, Clone, Copy)]
pub struct Line {
    pub weight: f64,
    pub slope: f64,
    pub reference: f64,
}

// Each contig's first half looks for its second half among every contig's second half. The weight
// that finds the most within NEIGHBOURS is the one that separates genomes on this assembly.
pub fn derive(contigs: &Contigs, presence_fraction: f64, seed: u64) -> Option<Derived> {
    let n = contigs.coverage.len();
    if n <= NEIGHBOURS {
        return None;
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let partners = noise::partners(&contigs.coverage, &contigs.whole, NEIGHBOURS, &mut rng)?;
    let found = (0..n)
        .into_par_iter()
        .map(|contig| recalled(contigs, &partners[contig], contig, presence_fraction))
        .collect::<Vec<_>>();
    let recall = (0..STEPS)
        .map(|step| found.iter().map(|held| held[step]).sum::<f64>() / n as f64)
        .collect::<Vec<_>>();
    let best = recall.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let plateau = (0..STEPS)
        .filter(|step| recall[*step] >= best - PLATEAU)
        .map(weight_at)
        .collect::<Vec<_>>();
    Some(Derived {
        weight: plateau.iter().sum::<f64>() / plateau.len() as f64,
        line: line(&found, &contigs.lengths),
    })
}

// A short contig's k-mers are a small sample of its genome, so the share coverage should carry
// can move with length. Fitted jointly over every contig rather than per band, which holds dozens.
fn line(found: &[[f64; STEPS]], lengths: &[usize]) -> Line {
    let logs = lengths
        .iter()
        .map(|length| ((*length / 2).max(1) as f64).log10())
        .collect::<Vec<_>>();
    let mut sorted = logs.clone();
    sorted.sort_unstable_by(f64::total_cmp);
    let reference = sorted[sorted.len() / 2];
    let scored = (0..SLOPES)
        .flat_map(|at| {
            let slope = SLOPE_REACH * (2.0 * at as f64 / (SLOPES - 1) as f64 - 1.0);
            (0..STEPS).map(move |step| (weight_at(step), slope))
        })
        .map(|(weight, slope)| {
            let recall = found
                .iter()
                .zip(&logs)
                .map(|(held, log)| held[nearest_step(weight + slope * (log - reference))])
                .sum::<f64>()
                / found.len() as f64;
            (recall, weight, slope)
        })
        .collect::<Vec<_>>();
    let best = scored
        .iter()
        .map(|(recall, _, _)| *recall)
        .fold(f64::NEG_INFINITY, f64::max);
    let plateau = scored
        .iter()
        .filter(|(recall, _, _)| *recall >= best - PLATEAU)
        .collect::<Vec<_>>();
    let mean = |pick: fn(&&(f64, f64, f64)) -> f64| {
        plateau.iter().map(pick).sum::<f64>() / plateau.len() as f64
    };
    Line {
        weight: mean(|point| point.1),
        slope: mean(|point| point.2),
        reference,
    }
}

fn nearest_step(weight: f64) -> usize {
    ((weight / STEP).round().max(0.0) as usize).min(STEPS - 1)
}

fn weight_at(step: usize) -> f64 {
    step as f64 * STEP
}

// The expectation over every candidate partner rather than one draw, which left the weight
// swinging twofold with the seed.
fn recalled(contigs: &Contigs, partners: &Partners, contig: usize, presence: f64) -> [f64; STEPS] {
    let coverage =
        |other: &[f64]| metabat_with(contigs.coverage[contig], other, MIN_VAR, MIN_VAR, presence).0;
    let own_composition = rho(&contigs.first[contig], &contigs.second[contig]);
    let own = partners
        .rows
        .iter()
        .map(|row| {
            let own_coverage = coverage(row);
            std::array::from_fn::<f64, STEPS, _>(|step| {
                combine(own_coverage, own_composition, weight_at(step))
            })
        })
        .collect::<Vec<_>>();
    let mut closer = vec![[0usize; STEPS]; own.len()];
    for other in (0..contigs.coverage.len()).filter(|other| *other != contig) {
        let distance = (
            coverage(contigs.coverage[other]),
            rho(&contigs.first[contig], &contigs.second[other]),
        );
        let at = std::array::from_fn::<f64, STEPS, _>(|step| {
            combine(distance.0, distance.1, weight_at(step))
        });
        for (counts, bars) in closer.iter_mut().zip(&own) {
            for step in 0..STEPS {
                counts[step] += usize::from(at[step] < bars[step]);
            }
        }
    }
    std::array::from_fn(|step| {
        closer
            .iter()
            .zip(&partners.weights)
            .map(|(counts, weight)| weight * f64::from(u8::from(counts[step] < NEIGHBOURS)))
            .sum()
    })
}
