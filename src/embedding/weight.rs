use rand::{SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::embedding::metrics::{DistanceSettings, combine, metabat_with, rho};

pub mod noise;

use noise::Partners;

pub const NEIGHBOURS: usize = 10;
const STEPS: usize = 20;
const STEP: f64 = 0.05;
const PLATEAU: f64 = 0.005;

pub struct Contigs<'a> {
    pub coverage: Vec<&'a [f64]>,
    pub whole: Vec<&'a [f64]>,
    pub first: Vec<Vec<f64>>,
    pub second: Vec<Vec<f64>>,
}

// Each contig's first half looks for its second half among every contig's second half. The weight
// that finds the most within NEIGHBOURS is the one that separates genomes on this assembly.
pub fn derive(contigs: &Contigs, distance: DistanceSettings, seed: u64) -> Option<f64> {
    let n = contigs.coverage.len();
    if n <= NEIGHBOURS {
        return None;
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let partners = noise::partners(&contigs.coverage, &contigs.whole, NEIGHBOURS, &mut rng)?;
    let found = (0..n)
        .into_par_iter()
        .map(|contig| recalled(contigs, &partners[contig], contig, distance))
        .collect::<Vec<_>>();
    let recall = (0..STEPS)
        .map(|step| found.iter().map(|held| held[step]).sum::<f64>() / n as f64)
        .collect::<Vec<_>>();
    let best = recall.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let plateau = (0..STEPS)
        .filter(|step| recall[*step] >= best - PLATEAU)
        .map(weight_at)
        .collect::<Vec<_>>();
    Some(plateau.iter().sum::<f64>() / plateau.len() as f64)
}

fn weight_at(step: usize) -> f64 {
    step as f64 * STEP
}

// The expectation over every candidate partner rather than one draw, which left the weight
// swinging twofold with the seed.
fn recalled(
    contigs: &Contigs,
    partners: &Partners,
    contig: usize,
    distance: DistanceSettings,
) -> [f64; STEPS] {
    let floor = distance.variance_floor;
    let coverage = |other: &[f64]| {
        metabat_with(
            contigs.coverage[contig],
            other,
            floor,
            floor,
            distance.presence_fraction,
        )
        .0
    };
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
