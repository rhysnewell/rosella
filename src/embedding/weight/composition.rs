use log::info;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::embedding::metrics::likelihood::{
    COMPOSITION_FLOOR, CompositionModel, Line, length_term, log_normal,
};
use crate::embedding::metrics::{Centred, rho_between};

use super::noise::{Density, RANDOM_PAIRS};

const MIN_HALVES: usize = 50;
const LOWEST_SPREAD: f64 = 0.01;
const STEPS: usize = 100;
const START_SHARE: f64 = 0.9;
const LOWEST_SHARE: f64 = 0.02;
// Median absolute deviation to the standard deviation of a normal.
const MAD_SCALE: f64 = 1.4826;

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

    let mut rng = StdRng::seed_from_u64(seed);
    let background = (0..RANDOM_PAIRS)
        .map(|_| (rng.random_range(0..n), rng.random_range(0..n)))
        .filter(|(a, b)| a != b)
        .map(|(a, b)| point(&centred[a], &centred[b], lengths[a], lengths[b]))
        .collect::<Vec<_>>();
    let other = least_squares(&background, &vec![1.0; background.len()])?;
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
            point(
                &Centred::new(halves.first[at]),
                &Centred::new(halves.second[at]),
                first,
                whole - first,
            )
        })
        .collect::<Vec<_>>();
    let same = least_squares(&same_points, &vec![1.0; same_points.len()])?;
    let spread = scatter(&same_points, &vec![1.0; same_points.len()], same);
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

pub fn point(a: &Centred, b: &Centred, a_length: f64, b_length: f64) -> (f64, f64) {
    (
        rho_between(a, b).max(COMPOSITION_FLOOR).ln(),
        length_term(a_length, b_length),
    )
}

// Pairs inside near complete bins are separate contigs of one genome, which halves are not, so
// they reset the same-genome line. The mixture against random pairs keeps contaminants out of it.
pub fn refit(model: &CompositionModel, points: &[(f64, f64)]) -> Option<CompositionModel> {
    if points.len() < MIN_HALVES {
        return None;
    }
    let other = points
        .iter()
        .map(|(x, u)| model.background.at(x - model.other.at(*u)).exp())
        .collect::<Vec<_>>();
    let mut weights = vec![1.0; points.len()];
    let mut same = least_squares(points, &weights)?;
    let mut residuals = points
        .iter()
        .map(|(x, u)| (x - same.at(*u)).abs())
        .collect::<Vec<_>>();
    let middle = residuals.len() / 2;
    let mut spread = (MAD_SCALE * *residuals.select_nth_unstable_by(middle, f64::total_cmp).1)
        .max(LOWEST_SPREAD);
    let mut share = START_SHARE;
    for _ in 0..STEPS {
        for ((weight, (x, u)), g) in weights.iter_mut().zip(points).zip(&other) {
            let f = log_normal(x - same.at(*u), spread * spread).exp();
            *weight = share * f / (share * f + (1.0 - share) * g);
        }
        share = (weights.iter().sum::<f64>() / points.len() as f64)
            .clamp(LOWEST_SHARE, 1.0 - LOWEST_SHARE);
        let Some(line) = least_squares(points, &weights) else {
            break;
        };
        same = line;
        spread = scatter(points, &weights, same);
    }
    info!(
        "Composition refit on {} pairs in near complete bins: share {share:.2}, same genome log \
         distance {:.2} + {:.2} x length term, spread {spread:.2}.",
        points.len(),
        same.intercept,
        same.slope,
    );
    Some(CompositionModel {
        same,
        spread,
        other: model.other,
        background: model.background.clone(),
    })
}

fn scatter(points: &[(f64, f64)], weights: &[f64], line: Line) -> f64 {
    let total = weights.iter().sum::<f64>();
    let squares = points
        .iter()
        .zip(weights)
        .map(|((x, u), weight)| weight * (x - line.at(*u)).powi(2))
        .sum::<f64>();
    (squares / total).sqrt().max(LOWEST_SPREAD)
}

fn least_squares(points: &[(f64, f64)], weights: &[f64]) -> Option<Line> {
    let total = weights.iter().sum::<f64>();
    if total <= 0.0 {
        return None;
    }
    let (mean_x, mean_u) = points
        .iter()
        .zip(weights)
        .fold((0.0, 0.0), |(x, u), ((px, pu), w)| (x + w * px, u + w * pu));
    let (mean_x, mean_u) = (mean_x / total, mean_u / total);
    let (covariance, variance) =
        points
            .iter()
            .zip(weights)
            .fold((0.0, 0.0), |(c, v), ((x, u), w)| {
                (
                    c + w * (u - mean_u) * (x - mean_x),
                    v + w * (u - mean_u).powi(2),
                )
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
