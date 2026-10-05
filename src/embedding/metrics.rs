use std::sync::Arc;

use erfc::erfc;
use likelihood::{Evidence, LogDepth};

pub mod erfc;
pub mod likelihood;
pub mod prepared;

const EPSILON: f64 = 1e-6;
pub use crate::tuning::MIN_VAR;
const MIN_VAR_EPSILON: f64 = 1e-4;
const SQRT_2: f64 = std::f64::consts::SQRT_2;

fn normal_cdf(mean: f64, sigma: f64, x: f64) -> f64 {
    (0.5 * erfc(-(x - mean) / (sigma * SQRT_2))).min(1.0)
}

// Folded rather than collected because this runs once per pairwise distance, which made
// the vector it replaces the program's hottest allocation.
#[derive(Default)]
struct Overlaps {
    total: f64,
    scored: usize,
}

impl Overlaps {
    fn push(&mut self, overlap: f64) {
        self.total += overlap;
        self.scored += 1;
    }

    // Nothing scored means both contigs are absent in every sample, which is agreement.
    fn finish(&self) -> (f64, usize) {
        if self.scored == 0 {
            return (EPSILON, 0);
        }
        let distance = self.total / self.scored as f64;
        (if distance.is_nan() { 1.0 } else { distance }, self.scored)
    }
}

// Both terms keep their own scale, so agreeing on one does not erase the other.
pub fn combine(coverage: f64, composition: f64, weight: f64) -> f64 {
    weight * coverage + (1.0 - weight) * composition
}

// The parts of the distance that are swept rather than derived. Carried as one value
// because the embedding and the refiner both compute it and must not drift apart.
#[derive(Debug, Clone, Default)]
pub struct DistanceSettings {
    pub presence_fraction: f64,
    pub aggregate_weight: Option<f64>,
    pub variance_unit: Option<f64>,
    pub evidence: Option<Arc<Evidence>>,
}

impl DistanceSettings {
    pub fn composition_only(&self) -> Self {
        Self {
            aggregate_weight: Some(0.0),
            evidence: None,
            ..self.clone()
        }
    }
}

// The mean shift, the variance clamp, the root and the log are all per row, so the prepared
// path hoists every one of them out of the pairwise loop.
#[derive(Debug, Clone, Copy)]
pub struct Moments {
    pub mean: f64,
    pub variance: f64,
    pub deviation: f64,
    pub log_variance: f64,
}

impl Moments {
    pub fn new(mean: f64, variance: f64) -> Self {
        Self {
            mean: mean + EPSILON,
            variance,
            deviation: variance.sqrt(),
            log_variance: variance.ln(),
        }
    }
}

pub fn overlap(a: Moments, b: Moments) -> f64 {
    // The neighbour join measures one pair from either end, and both have to round alike.
    let (a, b) = match (b.variance, b.mean) < (a.variance, a.mean) {
        true => (b, a),
        false => (a, b),
    };
    let (a_mean, a_var, a_sd) = (a.mean, a.variance, a.deviation);
    let (b_mean, b_var, b_sd) = (b.mean, b.variance, b.deviation);

    let (mut k1, mut k2) = if (a_var - b_var).abs() < MIN_VAR_EPSILON {
        let midpoint = (a_mean + b_mean) / 2.0;
        (midpoint, midpoint)
    } else {
        let tmp = (a_var
            * b_var
            * ((a_mean - b_mean) * (a_mean - b_mean)
                - (a_var - b_var) * (b.log_variance - a.log_variance)))
            .sqrt();
        (
            (tmp - a_mean * b_var + b_mean * a_var) / (a_var - b_var),
            (tmp + a_mean * b_var - b_mean * a_var) / (b_var - a_var),
        )
    };

    if k1 > k2 {
        std::mem::swap(&mut k1, &mut k2);
    }

    let ((narrow_mean, narrow_sd), (wide_mean, wide_sd)) = if a_var > b_var {
        ((b_mean, b_sd), (a_mean, a_sd))
    } else {
        ((a_mean, a_sd), (b_mean, b_sd))
    };

    if k1 == k2 {
        (normal_cdf(narrow_mean, narrow_sd, k1) - normal_cdf(wide_mean, wide_sd, k1)).abs()
    } else {
        (normal_cdf(narrow_mean, narrow_sd, k2) - normal_cdf(narrow_mean, narrow_sd, k1)
            + normal_cdf(wide_mean, wide_sd, k1)
            - normal_cdf(wide_mean, wide_sd, k2))
        .abs()
    }
}

// Mutual absence agrees for every pair of absent contigs, so scoring it would outvote the samples
// that saw something. Each contig's bar is a share of its own deepest sample.
pub(crate) fn coverage_distance(
    a: &[Moments],
    a_presence: f64,
    b: &[Moments],
    b_presence: f64,
) -> (f64, usize) {
    let mut overlaps = Overlaps::default();
    for (x, y) in a.iter().zip(b) {
        if x.mean - EPSILON <= a_presence && y.mean - EPSILON <= b_presence {
            continue;
        }
        overlaps.push(overlap(*x, *y).clamp(EPSILON, 1.0 - EPSILON));
    }
    overlaps.finish()
}

// Depth counted in units of `unit` is Poisson, so its variance is the depth itself, floored at one
// unit. The table's own column is per-base spread along the contig, not between contigs.
pub(crate) fn moments(row: &[f64], unit: Option<f64>) -> impl Iterator<Item = Moments> + '_ {
    row.chunks_exact(2).map(move |sample| {
        let variance = match unit {
            Some(unit) => unit * sample[0].max(unit),
            None => (sample[1] + EPSILON).max(MIN_VAR),
        };
        Moments::new(sample[0], variance)
    })
}

pub(crate) fn presence(row: &[f64], presence_fraction: f64) -> f64 {
    presence_fraction * peak_mean(row)
}

// What the abundance distance reads off one row, so a row met many times pays for it once.
pub struct Abundance {
    presence: f64,
    length: f64,
    samples: Samples,
}

enum Samples {
    Overlap(Vec<Moments>),
    Likelihood(Arc<Evidence>, Vec<LogDepth>),
}

impl Abundance {
    pub fn new(row: &[f64], length: usize, settings: &DistanceSettings) -> Self {
        let samples = match &settings.evidence {
            Some(evidence) => Samples::Likelihood(
                Arc::clone(evidence),
                evidence.depths(row, length as f64).collect(),
            ),
            None => Samples::Overlap(moments(row, settings.variance_unit).collect()),
        };
        Self {
            presence: presence(row, settings.presence_fraction),
            length: length as f64,
            samples,
        }
    }

    fn log_odds(&self, other: &Self) -> Option<(&Evidence, f64, usize)> {
        match (&self.samples, &other.samples) {
            (Samples::Likelihood(evidence, x), Samples::Likelihood(_, y)) => {
                let (odds, scored) = evidence.log_odds(x, self.presence, y, other.presence);
                Some((evidence, odds, scored))
            }
            _ => None,
        }
    }
}

pub fn abundance_distance(a: &Abundance, b: &Abundance) -> (f64, usize) {
    match (&a.samples, &b.samples) {
        (Samples::Overlap(x), Samples::Overlap(y)) => {
            coverage_distance(x, a.presence, y, b.presence)
        }
        _ => match a.log_odds(b) {
            Some((evidence, odds, scored)) => evidence.coverage(odds, scored),
            None => unreachable!("both rows are read under one distance settings"),
        },
    }
}

// The coverage distance alone and the distance both views give together.
pub fn aggregate(
    a: &Abundance,
    b: &Abundance,
    composition: f64,
    weight: Option<f64>,
) -> (f64, usize, f64) {
    if let Some((evidence, odds, scored)) = a.log_odds(b) {
        let together = odds + evidence.composition_ratio(composition, a.length, b.length);
        let (coverage, scored) = evidence.coverage(odds, scored);
        return (coverage, scored, evidence.posterior(together));
    }
    let (coverage, scored) = abundance_distance(a, b);
    let distance = combine(coverage, composition, weight_for(scored, weight));
    (
        coverage,
        scored,
        if distance.is_nan() { 1.0 } else { distance },
    )
}

fn peak_mean(row: &[f64]) -> f64 {
    row.iter()
        .step_by(2)
        .fold(0.0f64, |peak, mean| peak.max(*mean))
}

pub struct Centred {
    values: Vec<f64>,
    variance: f64,
}

impl Centred {
    pub fn new(row: &[f64]) -> Self {
        let (mean, variance) = centred_variance(row);
        Self {
            values: row.iter().map(|value| value - mean).collect(),
            variance,
        }
    }
}

// Proportionality distance. `vlr / (var(a) + var(b))`, which is `1 - rho`, on [0, 2].
pub fn rho_between(a: &Centred, b: &Centred) -> f64 {
    let covariance = a
        .values
        .iter()
        .zip(&b.values)
        .map(|(x, y)| x * y)
        .sum::<f64>();
    rho_from(covariance, a.variance, b.variance)
}

fn centred_variance(row: &[f64]) -> (f64, f64) {
    if row.is_empty() {
        return (0.0, 0.0);
    }
    let mean = row.iter().sum::<f64>() / row.len() as f64;
    let variance = row
        .iter()
        .map(|value| (value - mean) * (value - mean))
        .sum();
    (mean, variance)
}

pub(crate) fn rho_from(covariance: f64, var_a: f64, var_b: f64) -> f64 {
    let total_variance = var_a + var_b;
    if total_variance == 0.0 {
        return 0.0;
    }
    let distance = (-2.0 * covariance + total_variance) / total_variance;
    if distance.is_nan() { 2.0 } else { distance }
}

pub fn euclidean(a: &[f64], b: &[f64]) -> f64 {
    let distance = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt();
    if distance.is_nan() {
        f64::MAX
    } else {
        distance
    }
}

// Weight coverage against composition the way flight does, by sample count.
pub fn aggregate_weight(n_samples: usize) -> f64 {
    n_samples as f64 / (n_samples as f64 + 1.0)
}

pub fn weight_for(n_samples: usize, override_value: Option<f64>) -> f64 {
    override_value.unwrap_or_else(|| aggregate_weight(n_samples))
}

pub struct Point {
    pub abundance: Abundance,
    pub composition: Centred,
}

impl Point {
    pub fn new(
        coverage: &[f64],
        composition: &[f64],
        length: usize,
        settings: &DistanceSettings,
    ) -> Self {
        Self {
            abundance: Abundance::new(coverage, length, settings),
            composition: Centred::new(composition),
        }
    }
}

#[derive(Debug, Clone)]
pub struct AggregateMetric {
    n_coverage_columns: usize,
    settings: DistanceSettings,
}

impl AggregateMetric {
    pub fn new(n_coverage_columns: usize, settings: DistanceSettings) -> Self {
        Self {
            n_coverage_columns,
            settings,
        }
    }

    pub fn settings(&self) -> &DistanceSettings {
        &self.settings
    }

    pub fn distance(&self, a: &Point, b: &Point) -> f64 {
        aggregate(
            &a.abundance,
            &b.abundance,
            rho_between(&a.composition, &b.composition),
            self.settings.aggregate_weight,
        )
        .2
    }

    // Coverage is never negative and its weight never above the most the samples allow, so no
    // pair with this composition sits nearer. Rounding keeps the order, so the bound is exact.
    fn floor(&self, composition: f64) -> f64 {
        if composition <= 0.0 || composition.is_nan() {
            return f64::NEG_INFINITY;
        }
        let heaviest = weight_for(self.n_coverage_columns / 2, self.settings.aggregate_weight);
        combine(0.0, composition, heaviest)
    }

    fn combine(&self, coverage: f64, scored: usize, composition: f64) -> f64 {
        let weight = weight_for(scored, self.settings.aggregate_weight);
        let distance = combine(coverage, composition, weight);
        if distance.is_nan() { 1.0 } else { distance }
    }
}
