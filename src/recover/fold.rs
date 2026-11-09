use ndarray::{Array1, ArrayView2};
use rayon::prelude::*;

use crate::markers::replicon::DEPTH_FLOOR;
use crate::quality::{Bars, Quality};

const MARKER_FREE: f64 = 5.0;
// MIMAG medium quality (Bowers et al. 2017).
pub(crate) const HOST: Bars = Bars {
    completeness: 50.0,
    contamination: 10.0,
};
const CONFIDENT: f64 = 0.95;
pub(crate) const FITS: f64 = 0.05;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Role {
    Donor,
    Host,
    Other,
}

impl Role {
    pub fn of(quality: Quality, eukaryotic: bool) -> Self {
        if eukaryotic {
            Self::Other
        } else if quality.completeness < MARKER_FREE {
            Self::Donor
        } else if quality.clears(HOST) {
            Self::Host
        } else {
            Self::Other
        }
    }

    pub fn partial(quality: Quality) -> bool {
        quality.completeness < HOST.completeness
    }
}

struct Profile {
    mean: Array1<f64>,
    spread: Array1<f64>,
    weight: f64,
}

fn profile(members: &[usize], depths: ArrayView2<f64>, lengths: &[usize]) -> Profile {
    let logs = members
        .iter()
        .map(|contig| depths.row(*contig).mapv(|depth| (depth + DEPTH_FLOOR).ln()))
        .collect::<Vec<_>>();
    let weights = members
        .iter()
        .map(|contig| lengths[*contig] as f64)
        .collect::<Vec<_>>();
    let total = weights.iter().sum::<f64>();
    let mut mean = Array1::zeros(depths.ncols());
    for (log, weight) in logs.iter().zip(&weights) {
        mean.scaled_add(weight / total, log);
    }
    let mut spread = Array1::zeros(depths.ncols());
    for (log, weight) in logs.iter().zip(&weights) {
        spread.scaled_add(weight / total, &(log - &mean).mapv(|gap| gap * gap));
    }
    Profile {
        mean,
        spread,
        weight: total * total / weights.iter().map(|weight| weight * weight).sum::<f64>(),
    }
}

pub struct DepthFit {
    profiles: Vec<Profile>,
    shrunk: Vec<Array1<f64>>,
    samples: Vec<usize>,
}

impl DepthFit {
    pub fn of(
        groups: &[&[usize]],
        roles: &[Role],
        depths: ArrayView2<f64>,
        lengths: &[usize],
    ) -> Option<Self> {
        let hosts = (0..groups.len())
            .filter(|at| roles[*at] == Role::Host)
            .collect::<Vec<_>>();
        if hosts.is_empty() {
            return None;
        }
        let profiles = groups
            .par_iter()
            .map(|members| profile(members, depths, lengths))
            .collect::<Vec<_>>();
        let pooled = Array1::from_iter((0..depths.ncols()).map(|sample| {
            median(
                hosts
                    .iter()
                    .map(|host| profiles[*host].spread[sample])
                    .collect(),
            )
        }));
        let samples = pooled
            .iter()
            .enumerate()
            .filter(|(_, spread)| **spread > 0.0)
            .map(|(sample, _)| sample)
            .collect::<Vec<_>>();
        if samples.is_empty() {
            return None;
        }
        let shrunk = profiles
            .iter()
            .map(|profile| (&profile.spread * profile.weight + &pooled) / (profile.weight + 1.0))
            .collect();
        Some(Self {
            profiles,
            shrunk,
            samples,
        })
    }

    fn compare(&self, one: usize, other: usize) -> (f64, f64) {
        let (mut spread_sum, mut distance) = (0.0, 0.0);
        for sample in &self.samples {
            let spread = self.shrunk[other][*sample] + self.shrunk[one][*sample];
            let gap = self.profiles[one].mean[*sample] - self.profiles[other].mean[*sample];
            spread_sum += spread.ln();
            distance += gap * gap / spread;
        }
        (-0.5 * (spread_sum + distance), distance)
    }

    pub fn fits(&self, one: usize, other: usize) -> bool {
        chi_square_above(self.compare(one, other).1, self.samples.len()) >= FITS
    }

    // A donor reads as one contig of each group, and every group competes for it, so a sample many
    // genomes share at one depth spreads the posterior instead of handing it to the nearest.
    pub fn folds(&self, roles: &[Role]) -> Vec<(usize, usize)> {
        (0..roles.len())
            .into_par_iter()
            .filter(|donor| roles[*donor] == Role::Donor)
            .filter_map(|donor| {
                let fits = (0..roles.len())
                    .filter(|other| *other != donor)
                    .map(|other| (other, self.compare(donor, other).0))
                    .collect::<Vec<_>>();
                let (top, best) = fits
                    .iter()
                    .copied()
                    .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0)))?;
                let share = 1.0 / fits.iter().map(|(_, fit)| (fit - best).exp()).sum::<f64>();
                (roles[top] == Role::Host && share >= CONFIDENT && self.fits(donor, top))
                    .then_some((donor, top))
            })
            .collect()
    }
}

pub(crate) fn median(mut values: Vec<f64>) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_unstable_by(f64::total_cmp);
    let middle = values.len() / 2;
    match values.len() % 2 {
        1 => values[middle],
        _ => (values[middle - 1] + values[middle]) / 2.0,
    }
}

pub fn chi_square_above(statistic: f64, degrees: usize) -> f64 {
    let (a, x) = (degrees as f64 / 2.0, statistic / 2.0);
    if x <= 0.0 {
        return 1.0;
    }
    let front = (a * x.ln() - x - libm::lgamma(a)).exp();
    if x < a + 1.0 {
        let (mut term, mut sum, mut n) = (1.0 / a, 1.0 / a, a);
        while term > sum * f64::EPSILON {
            n += 1.0;
            term *= x / n;
            sum += term;
        }
        return (1.0 - sum * front).max(0.0);
    }
    let tiny = f64::MIN_POSITIVE / f64::EPSILON;
    let mut b = x + 1.0 - a;
    let (mut c, mut d) = (1.0 / tiny, 1.0 / b);
    let mut fraction = d;
    for step in 1..1000 {
        let an = -(step as f64) * (step as f64 - a);
        b += 2.0;
        d = an * d + b;
        d = if d.abs() < tiny { tiny } else { d };
        c = b + an / c;
        c = if c.abs() < tiny { tiny } else { c };
        d = 1.0 / d;
        fraction *= d * c;
        if (d * c - 1.0).abs() < f64::EPSILON {
            break;
        }
    }
    front * fraction
}
