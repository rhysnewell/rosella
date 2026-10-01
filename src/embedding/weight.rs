use rand::{SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::embedding::metrics::{Abundance, Centred, abundance_distance, combine, rho_between};

pub mod noise;

use noise::Partners;

pub const NEIGHBOURS: usize = 10;
pub const STEPS: usize = 20;
const STEP: f64 = 0.05;
const PLATEAU: f64 = 0.005;

pub struct Contigs<'a> {
    pub coverage: Vec<&'a [f64]>,
    pub whole: Vec<&'a [f64]>,
    pub first: Vec<&'a [f64]>,
    pub second: Vec<&'a [f64]>,
}

// Each contig's first half looks for its second half among every contig's second half. The weight
// that finds the most within NEIGHBOURS is the one that separates genomes on this assembly.
pub fn recall(contigs: &Contigs, presence_fraction: f64, seed: u64) -> Option<[f64; STEPS]> {
    let n = contigs.coverage.len();
    if n <= NEIGHBOURS {
        return None;
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let partners = noise::partners(&contigs.coverage, &contigs.whole, NEIGHBOURS, &mut rng)?;
    let prepared = Prepared {
        abundance: contigs
            .coverage
            .par_iter()
            .map(|row| Abundance::new(row, presence_fraction))
            .collect(),
        first: contigs
            .first
            .par_iter()
            .map(|row| Centred::new(row))
            .collect(),
        second: contigs
            .second
            .par_iter()
            .map(|row| Centred::new(row))
            .collect(),
        presence_fraction,
    };
    let found = (0..n)
        .into_par_iter()
        .map(|contig| recalled(&prepared, &partners[contig], contig))
        .collect::<Vec<_>>();
    Some(std::array::from_fn(|step| {
        found.iter().map(|held| held[step]).sum::<f64>() / n as f64
    }))
}

pub fn centre(recall: &[f64; STEPS]) -> f64 {
    let best = recall.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let plateau = (0..STEPS)
        .filter(|step| recall[*step] >= best - PLATEAU)
        .map(weight_at)
        .collect::<Vec<_>>();
    plateau.iter().sum::<f64>() / plateau.len() as f64
}

fn weight_at(step: usize) -> f64 {
    step as f64 * STEP
}

struct Prepared {
    abundance: Vec<Abundance>,
    first: Vec<Centred>,
    second: Vec<Centred>,
    presence_fraction: f64,
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
            let partner = Abundance::new(row, prepared.presence_fraction);
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
