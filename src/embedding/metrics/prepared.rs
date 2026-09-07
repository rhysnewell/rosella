use super::{
    AggregateMetric, CompositionMetric, CoverageAggregation, DistanceSettings, EPSILON, Moments,
    Overlaps, finish, overlap, peak_mean,
};
use crate::embedding::bands::DepthBands;

/// One flat buffer, and the composition half centred once as `f32`. The descent walks pairs in
/// no order it can prefetch, so a vector per row costs a cache miss on the only wide loop.
pub struct PreparedAggregate<'a> {
    metric: AggregateMetric<'a>,
    aggregation: CoverageAggregation,
    n_samples: usize,
    tnf_width: usize,
    samples: Vec<Moments>,
    presence: Vec<f64>,
    radii: Vec<f64>,
    keeps_weight: bool,
    composition: CompositionMetric,
    composition_scale: f64,
    tnf: Vec<f32>,
    tnf_variance: Vec<f32>,
}

impl<'a> PreparedAggregate<'a> {
    pub fn new(
        rows: &[Vec<f64>],
        floors: &[f64],
        n_coverage_columns: usize,
        settings: DistanceSettings,
        bands: Option<&'a DepthBands>,
    ) -> Self {
        let n_samples = n_coverage_columns / 2;
        let tnf_width = rows
            .first()
            .map_or(0, |row| row.len().saturating_sub(n_coverage_columns));

        let mut samples = Vec::with_capacity(rows.len() * n_samples);
        let mut presence = Vec::with_capacity(rows.len());
        let mut radii = Vec::with_capacity(match bands {
            Some(_) => rows.len() * n_samples,
            None => 0,
        });
        let mut tnf = Vec::with_capacity(rows.len() * tnf_width);
        let mut tnf_variance = Vec::with_capacity(rows.len());

        for (row, floor) in rows.iter().zip(floors) {
            let (coverage, composition) = row.split_at(n_coverage_columns);

            samples.extend(
                coverage
                    .chunks_exact(2)
                    .map(|sample| Moments::new(sample[0], (sample[1] + EPSILON).max(*floor))),
            );
            presence.push(settings.presence_fraction * peak_mean(coverage));
            if let Some(bands) = bands {
                radii.extend(
                    coverage
                        .chunks_exact(2)
                        .enumerate()
                        .map(|(sample, values)| bands.radius(sample, values[0])),
                );
            }

            let mean = if composition.is_empty() || !settings.composition.centres_rows() {
                0.0
            } else {
                composition.iter().sum::<f64>() / composition.len() as f64
            };
            let start = tnf.len();
            tnf.extend(composition.iter().map(|value| (value - mean) as f32));
            tnf_variance.push(dot(&tnf[start..], &tnf[start..]));
        }

        Self {
            metric: AggregateMetric::new(n_coverage_columns, settings).with_bands(bands),
            aggregation: settings.aggregation,
            n_samples,
            tnf_width,
            samples,
            presence,
            radii,
            keeps_weight: bands.is_some_and(DepthBands::keeps_weight),
            composition: settings.composition,
            composition_scale: settings.composition_scale,
            tnf,
            tnf_variance,
        }
    }

    pub fn distance(&self, a: usize, b: usize) -> f64 {
        let (coverage, scored) = self.coverage(a, b);
        self.metric
            .combine(coverage, scored, self.composition(a, b))
    }

    fn coverage(&self, a: usize, b: usize) -> (f64, usize) {
        let mut overlaps = Overlaps::new(self.aggregation);
        for (sample, (x, y)) in self
            .samples_of(a)
            .iter()
            .zip(self.samples_of(b))
            .enumerate()
        {
            if x.mean - EPSILON <= self.presence[a] && y.mean - EPSILON <= self.presence[b] {
                continue;
            }
            overlaps.seen += 1;
            if self.unresolved(a, b, sample, x.mean, y.mean) {
                continue;
            }
            overlaps.push(overlap(*x, *y).clamp(EPSILON, 1.0 - EPSILON));
        }
        finish(&overlaps, self.keeps_weight)
    }

    fn unresolved(&self, a: usize, b: usize, sample: usize, a_mean: f64, b_mean: f64) -> bool {
        if self.radii.is_empty() {
            return false;
        }
        let radius =
            self.radii[a * self.n_samples + sample].max(self.radii[b * self.n_samples + sample]);
        (a_mean - b_mean).abs() <= radius
    }

    fn composition(&self, a: usize, b: usize) -> f64 {
        self.composition.from_moments(
            dot(self.tnf_of(a), self.tnf_of(b)) as f64,
            self.tnf_variance[a] as f64,
            self.tnf_variance[b] as f64,
            self.composition_scale,
        )
    }

    fn samples_of(&self, row: usize) -> &[Moments] {
        &self.samples[row * self.n_samples..(row + 1) * self.n_samples]
    }

    fn tnf_of(&self, row: usize) -> &[f32] {
        &self.tnf[row * self.tnf_width..(row + 1) * self.tnf_width]
    }
}

fn dot(a: &[f32], b: &[f32]) -> f32 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
