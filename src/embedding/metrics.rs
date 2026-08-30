use itertools::izip;
use statrs::function::erf::erfc;

const EPSILON: f64 = 1e-6;
pub const MIN_VAR: f64 = 1.0;
const MIN_VAR_EPSILON: f64 = 1e-4;
const SQRT_2: f64 = std::f64::consts::SQRT_2;

/// How far the length-scaled variance floor is allowed to move from `MIN_VAR`. An unclamped
/// floor lets a megabase contig reach a variance near zero, which makes its coverage
/// distribution so sharp that everything else is maximally distant from it.
const VARIANCE_SCALE_RANGE: (f64, f64) = (0.25, 2.0);

fn normal_cdf(mean: f64, sigma: f64, x: f64) -> f64 {
    (0.5 * erfc(-(x - mean) / (sigma * SQRT_2))).min(1.0)
}

/// How the per-sample coverage distances become one number.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum CoverageAggregation {
    /// flight's. Dominated by its smallest term, so one agreeing sample pulls a pair
    /// together while every other sample disagrees.
    Geometric,
    #[default]
    Arithmetic,
    /// The worst sample decides, so a pair has to agree everywhere to be close.
    Max,
}

pub const AGGREGATION_NAMES: [&str; 3] = ["geometric", "arithmetic", "max"];

impl CoverageAggregation {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "geometric" => Some(Self::Geometric),
            "arithmetic" => Some(Self::Arithmetic),
            "max" => Some(Self::Max),
            _ => None,
        }
    }

    fn combine(&self, overlaps: &[f64]) -> f64 {
        match self {
            Self::Geometric => {
                (overlaps.iter().map(|d| d.ln()).sum::<f64>() / overlaps.len() as f64).exp()
            }
            Self::Arithmetic => overlaps.iter().sum::<f64>() / overlaps.len() as f64,
            Self::Max => overlaps.iter().copied().fold(f64::NAN, f64::max),
        }
    }
}

/// How coverage and composition become one number.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Combination {
    /// flight's. `metabat_with` clamps to `EPSILON`, so a coverage-agreeing pair lands three
    /// orders of magnitude below a composition-agreeing one and composition cannot outvote it.
    Geometric,
    /// Both terms keep their own scale, so agreeing on one does not erase the other.
    #[default]
    Arithmetic,
}

pub const COMBINATION_NAMES: [&str; 2] = ["geometric", "arithmetic"];

impl Combination {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "geometric" => Some(Self::Geometric),
            "arithmetic" => Some(Self::Arithmetic),
            _ => None,
        }
    }

    pub fn combine(&self, coverage: f64, composition: f64, weight: f64) -> f64 {
        match self {
            Self::Geometric => {
                (coverage.powf(weight) * composition.powf(1.0 - weight)).sqrt()
            }
            Self::Arithmetic => weight * coverage + (1.0 - weight) * composition,
        }
    }
}

/// The parts of the distance that are swept rather than derived. Carried as one value
/// because the embedding and the refiner both compute it and must not drift apart.
#[derive(Debug, Clone, Copy, Default)]
pub struct DistanceSettings {
    pub aggregation: CoverageAggregation,
    pub length_scaled_variance: bool,
    pub views: Views,
    pub aggregate_weight: Option<f64>,
    pub combination: Combination,
}

/// A contig's coverage is averaged over its own bases, so a long one is measured more
/// precisely. The flat `MIN_VAR` floor asserts the opposite for everything from 1.5 kb up.
pub fn variance_floor(length: usize, reference_length: usize, enabled: bool) -> f64 {
    if !enabled {
        return MIN_VAR;
    }
    let scale = (reference_length.max(1) as f64 / length.max(1) as f64).sqrt();
    MIN_VAR * scale.clamp(VARIANCE_SCALE_RANGE.0, VARIANCE_SCALE_RANGE.1)
}

/// MetaBAT abundance distance over a row of interleaved per-sample mean and variance.
///
/// flight skipped samples whose means agreed, so two contigs that agreed everywhere had
/// nothing left to average and came back maximally distant. Agreement is the strongest
/// evidence they share a genome, so those samples are scored like any other.
pub fn metabat_with(
    a: &[f64],
    b: &[f64],
    a_floor: f64,
    b_floor: f64,
    aggregation: CoverageAggregation,
) -> f64 {
    let n_samples = a.len() / 2;
    let mut overlaps = Vec::with_capacity(n_samples);

    let a_means = a.iter().step_by(2);
    let b_means = b.iter().step_by(2);
    let a_vars = a.iter().skip(1).step_by(2);
    let b_vars = b.iter().skip(1).step_by(2);

    for (a_mean, b_mean, a_var, b_var) in izip!(a_means, b_means, a_vars, b_vars) {
        let a_mean = a_mean + EPSILON;
        let b_mean = b_mean + EPSILON;
        let a_var = (a_var + EPSILON).max(a_floor);
        let b_var = (b_var + EPSILON).max(b_floor);

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

    let distance = aggregation.combine(&overlaps);
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

pub fn weight_for(n_samples: usize, override_value: Option<f64>) -> f64 {
    override_value.unwrap_or_else(|| aggregate_weight(n_samples))
}

/// Coverage and composition in one metric, over rows laid out as
/// `[interleaved mean/var .., clr tetranucleotide frequencies ..]`.
#[derive(Debug, Clone, Copy)]
pub struct AggregateMetric {
    n_coverage_columns: usize,
    weight: f64,
    settings: DistanceSettings,
}

impl AggregateMetric {
    pub fn new(n_coverage_columns: usize, settings: DistanceSettings) -> Self {
        Self {
            n_coverage_columns,
            weight: weight_for(n_coverage_columns / 2, settings.aggregate_weight),
            settings,
        }
    }

    pub fn distance(&self, a: &[f64], b: &[f64], a_floor: f64, b_floor: f64) -> f64 {
        let (a_coverage, a_tnf) = a.split_at(self.n_coverage_columns);
        let (b_coverage, b_tnf) = b.split_at(self.n_coverage_columns);

        let coverage_distance = metabat_with(
            a_coverage,
            b_coverage,
            a_floor,
            b_floor,
            self.settings.aggregation,
        );
        let composition_distance = rho(a_tnf, b_tnf);

        let distance =
            self.settings
                .combination
                .combine(coverage_distance, composition_distance, self.weight);
        if distance.is_nan() { 1.0 } else { distance }
    }
}

/// One of the three distances flight embedded separately before intersecting them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum View {
    Coverage,
    Rho,
    Euclidean,
}


pub const VIEW_NAMES: [&str; 4] = ["combined", "coverage", "rho", "euclidean"];

/// Which views get their own graph. All false is the single combined metric, which keeps the
/// default path unchanged.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Views {
    pub coverage: bool,
    pub rho: bool,
    pub euclidean: bool,
}

impl Views {
    pub fn parse(names: &[String]) -> Option<Self> {
        let mut views = Self::default();
        for name in names {
            match name.as_str() {
                "combined" => return (names.len() == 1).then_some(Self::default()),
                "coverage" => views.coverage = true,
                "rho" => views.rho = true,
                "euclidean" => views.euclidean = true,
                _ => return None,
            }
        }
        Some(views)
    }

    pub fn selected(&self) -> Vec<View> {
        let mut selected = Vec::with_capacity(3);
        if self.coverage {
            selected.push(View::Coverage);
        }
        if self.rho {
            selected.push(View::Rho);
        }
        if self.euclidean {
            selected.push(View::Euclidean);
        }
        selected
    }
}

/// One view's distance over the same concatenated row `AggregateMetric` splits, so the three
/// graphs are built from one copy of the features rather than three.
#[derive(Debug, Clone, Copy)]
pub struct ViewMetric {
    n_coverage_columns: usize,
    view: View,
    aggregation: CoverageAggregation,
}

impl ViewMetric {
    pub fn new(n_coverage_columns: usize, view: View, aggregation: CoverageAggregation) -> Self {
        Self {
            n_coverage_columns,
            view,
            aggregation,
        }
    }

    pub fn distance(&self, a: &[f64], b: &[f64], a_floor: f64, b_floor: f64) -> f64 {
        let (a_coverage, a_tnf) = a.split_at(self.n_coverage_columns);
        let (b_coverage, b_tnf) = b.split_at(self.n_coverage_columns);
        let distance = match self.view {
            View::Coverage => {
                metabat_with(a_coverage, b_coverage, a_floor, b_floor, self.aggregation)
            }
            View::Rho => rho(a_tnf, b_tnf),
            View::Euclidean => euclidean(a_tnf, b_tnf),
        };
        if distance.is_nan() { 1.0 } else { distance }
    }
}
