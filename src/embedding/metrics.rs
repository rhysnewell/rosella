use itertools::izip;
use statrs::function::erf::erfc;

const EPSILON: f64 = 1e-6;
const MIN_VAR: f64 = 1.0;
const MIN_VAR_EPSILON: f64 = 1e-4;
const SQRT_2: f64 = std::f64::consts::SQRT_2;

fn normal_cdf(mean: f64, sigma: f64, x: f64) -> f64 {
    (0.5 * erfc(-(x - mean) / (sigma * SQRT_2))).min(1.0)
}

/// MetaBAT abundance distance over a row of interleaved per-sample mean and variance.
/// Geometric mean of the per-sample overlap of two normal distributions.
///
/// flight skipped samples whose means agreed, so two contigs that agreed everywhere had
/// nothing left to average and came back maximally distant. Agreement is the strongest
/// evidence they share a genome, so those samples are scored like any other.
pub fn metabat(a: &[f64], b: &[f64]) -> f64 {
    let n_samples = a.len() / 2;
    let mut overlaps = Vec::with_capacity(n_samples);

    let a_means = a.iter().step_by(2);
    let b_means = b.iter().step_by(2);
    let a_vars = a.iter().skip(1).step_by(2);
    let b_vars = b.iter().skip(1).step_by(2);

    for (a_mean, b_mean, a_var, b_var) in izip!(a_means, b_means, a_vars, b_vars) {
        let a_mean = a_mean + EPSILON;
        let b_mean = b_mean + EPSILON;
        let a_var = (a_var + EPSILON).max(MIN_VAR);
        let b_var = (b_var + EPSILON).max(MIN_VAR);

        let (mut k1, mut k2) = if (a_var - b_var).abs() < MIN_VAR_EPSILON {
            let midpoint = (a_mean + b_mean) / 2.0;
            (midpoint, midpoint)
        } else {
            let tmp = (a_var
                * b_var
                * ((a_mean - b_mean) * (a_mean - b_mean)
                    - 2.0 * (a_var - b_var) * (b_var / a_var).sqrt().ln()))
            .sqrt();
            (
                (tmp - a_mean * b_var + b_mean * a_var) / (a_var - b_var),
                (tmp + a_mean * b_var - b_mean * a_var) / (b_var - a_var),
            )
        };

        if k1 > k2 {
            std::mem::swap(&mut k1, &mut k2);
        }

        let ((narrow_mean, narrow_var), (wide_mean, wide_var)) = if a_var > b_var {
            ((b_mean, b_var), (a_mean, a_var))
        } else {
            ((a_mean, a_var), (b_mean, b_var))
        };
        let narrow_sd = narrow_var.sqrt();
        let wide_sd = wide_var.sqrt();

        let overlap = if k1 == k2 {
            (normal_cdf(narrow_mean, narrow_sd, k1) - normal_cdf(wide_mean, wide_sd, k1)).abs()
        } else {
            (normal_cdf(narrow_mean, narrow_sd, k2) - normal_cdf(narrow_mean, narrow_sd, k1)
                + normal_cdf(wide_mean, wide_sd, k1)
                - normal_cdf(wide_mean, wide_sd, k2))
            .abs()
        };

        // An unclamped zero drags the geometric mean to zero on its own.
        overlaps.push(overlap.clamp(EPSILON, 1.0 - EPSILON));
    }

    if overlaps.is_empty() {
        return 1.0;
    }

    let log_mean = overlaps.iter().map(|d| d.ln()).sum::<f64>() / overlaps.len() as f64;
    let distance = log_mean.exp();
    if distance.is_nan() { 1.0 } else { distance }
}

/// Proportionality distance. `vlr / (var(a) + var(b))`, which is `1 - rho`, on [0, 2].
pub fn rho(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let mean_a = a.iter().sum::<f64>() / n;
    let mean_b = b.iter().sum::<f64>() / n;

    let mut var_a = 0.0;
    let mut var_b = 0.0;
    let mut covariance = 0.0;
    for (x, y) in a.iter().zip(b.iter()) {
        let x = x - mean_a;
        let y = y - mean_b;
        var_a += x * x;
        var_b += y * y;
        covariance += x * y;
    }

    let total_variance = var_a + var_b;
    if total_variance == 0.0 {
        return 0.0;
    }

    let log_ratio_variance = -2.0 * covariance + var_a + var_b;
    let distance = log_ratio_variance / total_variance;
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

/// Weight coverage against composition the way flight does, by sample count.
pub fn aggregate_weight(n_samples: usize) -> f64 {
    n_samples as f64 / (n_samples as f64 + 1.0)
}

/// Coverage and composition in one metric, over rows laid out as
/// `[interleaved mean/var .., clr tetranucleotide frequencies ..]`.
#[derive(Debug, Clone, Copy)]
pub struct AggregateMetric {
    n_coverage_columns: usize,
    weight: f64,
}

impl AggregateMetric {
    pub fn new(n_coverage_columns: usize) -> Self {
        Self {
            n_coverage_columns,
            weight: aggregate_weight(n_coverage_columns / 2),
        }
    }

    pub fn distance(&self, a: &[f64], b: &[f64]) -> f64 {
        let (a_coverage, a_tnf) = a.split_at(self.n_coverage_columns);
        let (b_coverage, b_tnf) = b.split_at(self.n_coverage_columns);

        let coverage_distance = metabat(a_coverage, b_coverage);
        let composition_distance = rho(a_tnf, b_tnf);

        let distance = (coverage_distance.powf(self.weight)
            * composition_distance.powf(1.0 - self.weight))
        .sqrt();
        if distance.is_nan() { 1.0 } else { distance }
    }
}
