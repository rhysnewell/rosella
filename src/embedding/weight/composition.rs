use log::info;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::embedding::metrics::likelihood::{
    COMPOSITION_FLOOR, CompositionModel, Line, TAILS, Tail, length_term,
};
use crate::embedding::metrics::{Centred, rho_between};

use super::noise::{Density, RANDOM_PAIRS};

const MIN_HALVES: usize = 50;
const LOWEST_SCALE: f64 = 1e-4;
const SCALE_STEPS: usize = 100;

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
    let other = least_squares(&background)?;
    let density = Density::fit(
        &background
            .iter()
            .map(|(x, u)| x - other.at(*u))
            .collect::<Vec<_>>(),
    )?;

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
    let same = least_squares(&same_points)?;
    let residuals = same_points
        .iter()
        .map(|(x, u)| x - same.at(*u))
        .collect::<Vec<_>>();
    let (_, scale, tail) = TAILS
        .iter()
        .map(|degrees| student(&residuals, Tail::new(*degrees)))
        .max_by(|a, b| a.0.total_cmp(&b.0))?;
    info!(
        "Composition model from {} split contigs: same genome log distance {:.2} + {:.2} x length \
         term, scale {:.2}, tail {}, against random pairs {:.2} + {:.2} x length term.",
        same_points.len(),
        same.intercept,
        same.slope,
        scale.sqrt(),
        tail.degrees,
        other.intercept,
        other.slope
    );
    Some(CompositionModel {
        same,
        scale,
        tail,
        other,
        background: density.logged(),
    })
}

fn point(a: &Centred, b: &Centred, a_length: f64, b_length: f64) -> (f64, f64) {
    (
        rho_between(a, b).max(COMPOSITION_FLOOR).ln(),
        length_term(a_length, b_length),
    )
}

// The scale at a fixed tail by reweighting (Lange et al. 1989), each point down-weighted by how
// far out in the tail it sits.
fn student(residuals: &[f64], tail: Tail) -> (f64, f64, Tail) {
    let n = residuals.len() as f64;
    let mut scale = (residuals.iter().map(|r| r * r).sum::<f64>() / n).max(LOWEST_SCALE);
    if tail.degrees.is_finite() {
        for _ in 0..SCALE_STEPS {
            scale = (residuals
                .iter()
                .map(|r| r * r * (tail.degrees + 1.0) / (tail.degrees + r * r / scale))
                .sum::<f64>()
                / n)
                .max(LOWEST_SCALE);
        }
    }
    let log_likelihood = residuals.iter().map(|r| tail.log_density(*r, scale)).sum();
    (log_likelihood, scale, tail)
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
