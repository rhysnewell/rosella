use log::info;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::embedding::metrics::likelihood::{
    COMPOSITION_FLOOR, CompositionModel, Line, length_term,
};
use crate::embedding::metrics::{Centred, rho_between};

use super::noise::{Density, RANDOM_PAIRS};

const MIN_HALVES: usize = 50;
const LOWEST_SPREAD: f64 = 0.01;

pub struct Halves<'a> {
    pub first: Vec<&'a [f64]>,
    pub second: Vec<&'a [f64]>,
    pub lengths: Vec<f64>,
}

// Two halves of one contig are the one sample of genome mates that needs no labels and reads no
// depth, so they set the same-genome line and random pairs set the other.
pub fn fit(
    lengths: &[f64],
    composition: &[&[f64]],
    halves: &Halves,
    seed: u64,
) -> Option<CompositionModel> {
    if halves.lengths.len() < MIN_HALVES {
        return None;
    }
    let n = composition.len();
    let centred = composition
        .par_iter()
        .map(|row| Centred::new(row))
        .collect::<Vec<_>>();
    let log_distance = |a: &Centred, b: &Centred| rho_between(a, b).max(COMPOSITION_FLOOR).ln();

    let mut rng = StdRng::seed_from_u64(seed);
    let background = (0..RANDOM_PAIRS)
        .map(|_| (rng.random_range(0..n), rng.random_range(0..n)))
        .filter(|(a, b)| a != b)
        .map(|(a, b)| {
            (
                log_distance(&centred[a], &centred[b]),
                length_term(lengths[a], lengths[b]),
            )
        })
        .collect::<Vec<_>>();
    let other = least_squares(&background)?;
    let residuals = background
        .iter()
        .map(|(x, u)| x - other.at(*u))
        .collect::<Vec<_>>();
    let density = Density::fit(&residuals)?;

    let same_points = (0..halves.lengths.len())
        .into_par_iter()
        .map(|at| {
            let whole = halves.lengths[at];
            let first = (whole / 2.0).floor();
            (
                log_distance(
                    &Centred::new(halves.first[at]),
                    &Centred::new(halves.second[at]),
                ),
                length_term(first, whole - first),
            )
        })
        .collect::<Vec<_>>();
    let same = least_squares(&same_points)?;
    let spread = (same_points
        .iter()
        .map(|(x, u)| (x - same.at(*u)).powi(2))
        .sum::<f64>()
        / same_points.len() as f64)
        .sqrt()
        .max(LOWEST_SPREAD);
    info!(
        "Composition model from {} split contigs: same genome log distance {:.2} + {:.2} x length \
         term, spread {spread:.2}, against random pairs {:.2} + {:.2} x length term.",
        same_points.len(),
        same.intercept,
        same.slope,
        other.intercept,
        other.slope
    );
    Some(CompositionModel {
        same,
        spread,
        other,
        background: density.logged(),
    })
}

fn least_squares(points: &[(f64, f64)]) -> Option<Line> {
    let n = points.len() as f64;
    let (mean_x, mean_u) = points.iter().fold((0.0, 0.0), |(x, u), point| {
        (x + point.0 / n, u + point.1 / n)
    });
    let (covariance, variance) = points.iter().fold((0.0, 0.0), |(c, v), (x, u)| {
        (c + (u - mean_u) * (x - mean_x), v + (u - mean_u).powi(2))
    });
    if variance <= f64::EPSILON {
        return None;
    }
    let slope = covariance / variance;
    Some(Line {
        intercept: mean_x - slope * mean_u,
        slope,
    })
}
