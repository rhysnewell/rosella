use rand::{SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::embedding::metrics::{
    Abundance, Centred, DistanceSettings, abundance_distance, combine, rho_between,
};

pub mod noise;

use noise::Partners;

pub const NEIGHBOURS: usize = 10;
pub const STEPS: usize = 20;
const STEP: f64 = 0.05;
const PLATEAU: f64 = 0.005;
pub const UNITS: usize = 13;
const LOWEST_UNIT: f64 = -2.0;
const UNIT_STEP: f64 = 0.25;

pub struct Contigs<'a> {
    pub coverage: Vec<&'a [f64]>,
    pub whole: Vec<&'a [f64]>,
    pub first: Vec<&'a [f64]>,
    pub second: Vec<&'a [f64]>,
}

// Each contig's first half looks for its second half among every contig's second half. The weight
// that finds the most within NEIGHBOURS is the one that separates genomes on this assembly.
pub fn recall(contigs: &Contigs, settings: DistanceSettings, seed: u64) -> Option<[f64; STEPS]> {
    recall_each(contigs, &[settings], seed).map(|curves| curves[0])
}

// The partners and the composition halves do not depend on the coverage distance, so every setting
// is judged against the same draw.
pub fn recall_each(
    contigs: &Contigs,
    settings: &[DistanceSettings],
    seed: u64,
) -> Option<Vec<[f64; STEPS]>> {
    let n = contigs.coverage.len();
    if n <= NEIGHBOURS {
        return None;
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let partners = noise::partners(&contigs.coverage, &contigs.whole, NEIGHBOURS, &mut rng)?;
    let halves = |rows: &[&[f64]]| -> Vec<Centred> {
        rows.par_iter().map(|row| Centred::new(row)).collect()
    };
    let (first, second) = (halves(&contigs.first), halves(&contigs.second));
    let curves = settings
        .iter()
        .map(|settings| {
            let prepared = Prepared {
                abundance: contigs
                    .coverage
                    .par_iter()
                    .map(|row| Abundance::new(row, *settings))
                    .collect(),
                first: &first,
                second: &second,
                settings: *settings,
            };
            let found = (0..n)
                .into_par_iter()
                .map(|contig| recalled(&prepared, &partners[contig], contig))
                .collect::<Vec<_>>();
            std::array::from_fn(|step| found.iter().map(|held| held[step]).sum::<f64>() / n as f64)
        })
        .collect();
    Some(curves)
}

pub fn centre(recall: &[f64; STEPS]) -> f64 {
    plateau_centre((0..STEPS).map(|step| (weight_at(step), recall[step])))
}

// A unit is judged at its own best weight, and the centre is taken on a log scale because the
// unit is a depth.
pub fn unit_centre(curves: &[[f64; STEPS]]) -> f64 {
    let best = |curve: &[f64; STEPS]| curve.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    10f64.powf(plateau_centre(
        curves
            .iter()
            .enumerate()
            .map(|(at, curve)| (unit_at(at).log10(), best(curve))),
    ))
}

fn plateau_centre(points: impl Iterator<Item = (f64, f64)>) -> f64 {
    let points = points.collect::<Vec<_>>();
    let best = points
        .iter()
        .map(|(_, recall)| *recall)
        .fold(f64::NEG_INFINITY, f64::max);
    let plateau = points
        .iter()
        .filter(|(_, recall)| *recall >= best - PLATEAU)
        .map(|(at, _)| *at)
        .collect::<Vec<_>>();
    plateau.iter().sum::<f64>() / plateau.len() as f64
}

fn weight_at(step: usize) -> f64 {
    step as f64 * STEP
}

// One read on a contig at the cutoff is about 0.1 depth for short reads and several for long
// reads, so the search spans both.
pub fn unit_at(step: usize) -> f64 {
    10f64.powf(LOWEST_UNIT + step as f64 * UNIT_STEP)
}

struct Prepared<'a> {
    abundance: Vec<Abundance>,
    first: &'a [Centred],
    second: &'a [Centred],
    settings: DistanceSettings,
}

// The expectation over every candidate partner rather than one draw, which left the weight
// swinging twofold with the seed.
fn recalled(prepared: &Prepared, partners: &Partners, contig: usize) -> [f64; STEPS] {
    let mine = &prepared.abundance[contig];
    let own_composition = rho_between(&prepared.first[contig], &prepared.second[contig]);
    let own = partners
        .rows
        .iter()
        .map(|row| {
            let partner = Abundance::new(row, prepared.settings);
            let own_coverage = abundance_distance(mine, &partner).0;
            std::array::from_fn::<f64, STEPS, _>(|step| {
                combine(own_coverage, own_composition, weight_at(step))
            })
        })
        .collect::<Vec<_>>();
    // A partner is recalled when fewer than NEIGHBOURS others sit nearer, which is whether the
    // NEIGHBOURS-th nearest other sits no nearer than it, so one selection answers every partner.
    let others = prepared.abundance.len() - 1;
    let mut columns: [Vec<f64>; STEPS] = std::array::from_fn(|_| Vec::with_capacity(others));
    for other in (0..prepared.abundance.len()).filter(|other| *other != contig) {
        let coverage = abundance_distance(mine, &prepared.abundance[other]).0;
        let composition = rho_between(&prepared.first[contig], &prepared.second[other]);
        for (step, column) in columns.iter_mut().enumerate() {
            column.push(combine(coverage, composition, weight_at(step)));
        }
    }
    let nearest = columns
        .iter_mut()
        .map(|column| {
            *column
                .select_nth_unstable_by(NEIGHBOURS - 1, f64::total_cmp)
                .1
        })
        .collect::<Vec<_>>();
    std::array::from_fn(|step| {
        own.iter()
            .zip(&partners.weights)
            .map(|(bars, weight)| weight * f64::from(u8::from(nearest[step] >= bars[step])))
            .sum()
    })
}
