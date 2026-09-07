use itertools::izip;

use crate::embedding::bands::DepthBands;
use statrs::function::erf::erfc;

pub mod prepared;

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
}

/// Folded rather than collected because this runs once per pairwise distance, which made
/// the vector it replaces the program's hottest allocation.
struct Overlaps {
    aggregation: CoverageAggregation,
    total: f64,
    scored: usize,
    seen: usize,
}

impl Overlaps {
    fn new(aggregation: CoverageAggregation) -> Self {
        let total = match aggregation {
            CoverageAggregation::Max => f64::NAN,
            _ => 0.0,
        };
        Self {
            aggregation,
            total,
            scored: 0,
            seen: 0,
        }
    }

    fn push(&mut self, overlap: f64) {
        self.total = match self.aggregation {
            CoverageAggregation::Geometric => self.total + overlap.ln(),
            CoverageAggregation::Arithmetic => self.total + overlap,
            CoverageAggregation::Max => self.total.max(overlap),
        };
        self.scored += 1;
    }

    fn finish(&self) -> f64 {
        match self.aggregation {
            CoverageAggregation::Geometric => (self.total / self.scored as f64).exp(),
            CoverageAggregation::Arithmetic => self.total / self.scored as f64,
            CoverageAggregation::Max => self.total,
        }
    }
}

/// What a sample that cannot resolve a pair does to the coverage half. `Drop` lets the pair
/// lose the sample's weight as well as its vote, `Keep` charges the vote only.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum CoverageBand {
    #[default]
    Off,
    Drop,
    Keep,
}

pub const BAND_NAMES: [&str; 3] = ["off", "drop", "keep"];

impl CoverageBand {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "off" => Some(Self::Off),
            "drop" => Some(Self::Drop),
            "keep" => Some(Self::Keep),
            _ => None,
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
            Self::Geometric => (coverage.powf(weight) * composition.powf(1.0 - weight)).sqrt(),
            Self::Arithmetic => weight * coverage + (1.0 - weight) * composition,
        }
    }
}

/// How two composition rows become one number. Every variant returns a distance on [0, 2],
/// because the refiner's rho and aggregate bars are calibrated to that range.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum CompositionMetric {
    #[default]
    Rho,
    Cosine,
    Aitchison,
    Hellinger,
    TetraZ,
}

pub const COMPOSITION_NAMES: [&str; 5] = ["rho", "cosine", "aitchison", "hellinger", "tetra-z"];

/// Hellinger is `||sqrt(p) - sqrt(q)|| / sqrt(2)` on [0, 1], and the bars want [0, 2].
const HELLINGER_SCALE: f64 = std::f64::consts::FRAC_1_SQRT_2;

impl CompositionMetric {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "rho" => Some(Self::Rho),
            "cosine" => Some(Self::Cosine),
            "aitchison" => Some(Self::Aitchison),
            "hellinger" => Some(Self::Hellinger),
            "tetra-z" => Some(Self::TetraZ),
            _ => None,
        }
    }

    /// The correlation family reads centred moments. The two straight distances must not be
    /// shifted at all, because moving each row to its own mean moves the gap between them.
    pub fn centres_rows(&self) -> bool {
        matches!(self, Self::Rho | Self::Cosine | Self::TetraZ)
    }

    /// Aitchison is the only variant whose spread depends on the table, so it is the only one
    /// that reads the run-derived scale.
    pub fn scale(&self, run_scale: f64) -> f64 {
        match self {
            Self::Aitchison if run_scale > 0.0 => run_scale,
            Self::Hellinger => HELLINGER_SCALE,
            _ => 1.0,
        }
    }

    pub fn distance(&self, a: &[f64], b: &[f64], run_scale: f64) -> f64 {
        match self {
            Self::Rho => rho(a, b),
            Self::Cosine | Self::TetraZ => correlation(a, b),
            Self::Aitchison | Self::Hellinger => scaled_l2(a, b, self.scale(run_scale)),
        }
    }

    pub fn from_moments(&self, covariance: f64, var_a: f64, var_b: f64, run_scale: f64) -> f64 {
        match self {
            Self::Rho => rho_from(covariance, var_a, var_b),
            Self::Cosine | Self::TetraZ => correlation_from(covariance, var_a, var_b),
            Self::Aitchison | Self::Hellinger => {
                scaled_l2_from(var_a + var_b - 2.0 * covariance, self.scale(run_scale))
            }
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
    pub presence_fraction: f64,
    pub coverage_band: CoverageBand,
    pub composition: CompositionMetric,
    pub composition_scale: f64,
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

/// The mean shift, the variance clamp, the root and the log are all per row, so the prepared
/// path hoists every one of them out of the pairwise loop.
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

fn overlap(a: Moments, b: Moments) -> f64 {
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

fn finish(overlaps: &Overlaps, weigh_by_seen: bool) -> (f64, usize) {
    // Nothing scored means both contigs are absent in every sample, which is agreement.
    if overlaps.scored == 0 {
        return (EPSILON, 0);
    }
    let counted = match weigh_by_seen {
        true => overlaps.seen,
        false => overlaps.scored,
    };
    let distance = overlaps.finish();
    (if distance.is_nan() { 1.0 } else { distance }, counted)
}

/// MetaBAT abundance distance, with the count of samples that carried evidence.
///
/// A sample where both contigs are absent agrees for every pair of absent contigs, so scoring it
/// lets mutual absence outvote the samples that saw something. flight skipped samples whose means
/// *agreed*, which is the opposite condition and throws away real evidence; those stay scored.
/// Each contig's bar is a fraction of its own deepest sample, so a dense table skips nothing and
/// a deep contig cannot mask a shallow partner that is genuinely there.
pub fn metabat_with(
    a: &[f64],
    b: &[f64],
    a_floor: f64,
    b_floor: f64,
    aggregation: CoverageAggregation,
    presence_fraction: f64,
    bands: Option<&DepthBands>,
) -> (f64, usize) {
    let mut overlaps = Overlaps::new(aggregation);

    let a_presence = presence_fraction * peak_mean(a);
    let b_presence = presence_fraction * peak_mean(b);

    let a_means = a.iter().step_by(2);
    let b_means = b.iter().step_by(2);
    let a_vars = a.iter().skip(1).step_by(2);
    let b_vars = b.iter().skip(1).step_by(2);

    for (sample, (a_mean, b_mean, a_var, b_var)) in
        izip!(a_means, b_means, a_vars, b_vars).enumerate()
    {
        if *a_mean <= a_presence && *b_mean <= b_presence {
            continue;
        }
        overlaps.seen += 1;
        if let Some(bands) = bands
            && unresolved(bands, sample, *a_mean, *b_mean)
        {
            continue;
        }
        let a_var = (a_var + EPSILON).max(a_floor);
        let b_var = (b_var + EPSILON).max(b_floor);

        // An unclamped zero drags the geometric mean to zero on its own.
        overlaps.push(
            overlap(Moments::new(*a_mean, a_var), Moments::new(*b_mean, b_var))
                .clamp(EPSILON, 1.0 - EPSILON),
        );
    }

    finish(&overlaps, bands.is_some_and(DepthBands::keeps_weight))
}

fn unresolved(bands: &DepthBands, sample: usize, a_mean: f64, b_mean: f64) -> bool {
    let radius = bands
        .radius(sample, a_mean)
        .max(bands.radius(sample, b_mean));
    (a_mean - b_mean).abs() <= radius
}

fn peak_mean(row: &[f64]) -> f64 {
    row.iter()
        .step_by(2)
        .fold(0.0f64, |peak, mean| peak.max(*mean))
}

/// Proportionality distance. `vlr / (var(a) + var(b))`, which is `1 - rho`, on [0, 2].
pub fn rho(a: &[f64], b: &[f64]) -> f64 {
    let (mean_a, var_a) = centred_variance(a);
    let (mean_b, var_b) = centred_variance(b);
    let covariance = a
        .iter()
        .zip(b)
        .map(|(x, y)| (x - mean_a) * (y - mean_b))
        .sum::<f64>();
    rho_from(covariance, var_a, var_b)
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

fn rho_from(covariance: f64, var_a: f64, var_b: f64) -> f64 {
    let total_variance = var_a + var_b;
    if total_variance == 0.0 {
        return 0.0;
    }
    let distance = (-2.0 * covariance + total_variance) / total_variance;
    if distance.is_nan() { 2.0 } else { distance }
}

/// Correlation distance, `1 - r`. rho's sibling: the same covariance over a geometric mean of
/// the two variances rather than an arithmetic one.
pub fn correlation(a: &[f64], b: &[f64]) -> f64 {
    let (mean_a, var_a) = centred_variance(a);
    let (mean_b, var_b) = centred_variance(b);
    let covariance = a
        .iter()
        .zip(b)
        .map(|(x, y)| (x - mean_a) * (y - mean_b))
        .sum::<f64>();
    correlation_from(covariance, var_a, var_b)
}

fn correlation_from(covariance: f64, var_a: f64, var_b: f64) -> f64 {
    let spread = (var_a * var_b).sqrt();
    if spread <= 0.0 {
        return 0.0;
    }
    let distance = 1.0 - covariance / spread;
    if distance.is_nan() {
        2.0
    } else {
        distance.clamp(0.0, 2.0)
    }
}

pub fn scaled_l2(a: &[f64], b: &[f64], scale: f64) -> f64 {
    let squared = a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum::<f64>();
    scaled_l2_from(squared, scale)
}

fn scaled_l2_from(squared: f64, scale: f64) -> f64 {
    let distance = squared.max(0.0).sqrt() / scale;
    if distance.is_nan() {
        2.0
    } else {
        distance.min(2.0)
    }
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
pub struct AggregateMetric<'a> {
    n_coverage_columns: usize,
    settings: DistanceSettings,
    bands: Option<&'a DepthBands>,
}

impl<'a> AggregateMetric<'a> {
    pub fn new(n_coverage_columns: usize, settings: DistanceSettings) -> Self {
        Self {
            n_coverage_columns,
            settings,
            bands: None,
        }
    }

    pub fn with_bands(mut self, bands: Option<&'a DepthBands>) -> Self {
        self.bands = bands;
        self
    }

    pub fn distance(&self, a: &[f64], b: &[f64], a_floor: f64, b_floor: f64) -> f64 {
        let (a_coverage, a_tnf) = a.split_at(self.n_coverage_columns);
        let (b_coverage, b_tnf) = b.split_at(self.n_coverage_columns);

        let (coverage_distance, scored) = metabat_with(
            a_coverage,
            b_coverage,
            a_floor,
            b_floor,
            self.settings.aggregation,
            self.settings.presence_fraction,
            self.bands,
        );
        let composition =
            self.settings
                .composition
                .distance(a_tnf, b_tnf, self.settings.composition_scale);
        self.combine(coverage_distance, scored, composition)
    }

    fn combine(&self, coverage: f64, scored: usize, composition: f64) -> f64 {
        let weight = weight_for(scored, self.settings.aggregate_weight);
        let distance = self
            .settings
            .combination
            .combine(coverage, composition, weight);
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

impl View {
    pub fn name(&self) -> &'static str {
        match self {
            Self::Coverage => VIEW_NAMES[1],
            Self::Rho => VIEW_NAMES[2],
            Self::Euclidean => VIEW_NAMES[3],
        }
    }
}

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
pub struct ViewMetric<'a> {
    n_coverage_columns: usize,
    view: View,
    settings: DistanceSettings,
    bands: Option<&'a DepthBands>,
}

impl<'a> ViewMetric<'a> {
    pub fn new(n_coverage_columns: usize, view: View, settings: DistanceSettings) -> Self {
        Self {
            n_coverage_columns,
            view,
            settings,
            bands: None,
        }
    }

    pub fn with_bands(mut self, bands: Option<&'a DepthBands>) -> Self {
        self.bands = bands;
        self
    }

    pub fn distance(&self, a: &[f64], b: &[f64], a_floor: f64, b_floor: f64) -> f64 {
        let (a_coverage, a_tnf) = a.split_at(self.n_coverage_columns);
        let (b_coverage, b_tnf) = b.split_at(self.n_coverage_columns);
        let distance = match self.view {
            View::Coverage => {
                metabat_with(
                    a_coverage,
                    b_coverage,
                    a_floor,
                    b_floor,
                    self.settings.aggregation,
                    self.settings.presence_fraction,
                    self.bands,
                )
                .0
            }
            View::Rho => {
                self.settings
                    .composition
                    .distance(a_tnf, b_tnf, self.settings.composition_scale)
            }
            View::Euclidean => euclidean(a_tnf, b_tnf),
        };
        if distance.is_nan() { 1.0 } else { distance }
    }
}
