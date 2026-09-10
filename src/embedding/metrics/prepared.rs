use ndarray::Array2;

use super::{
    AggregateMetric, Combination, CompositionMetric, CoverageAggregation, DistanceSettings,
    EPSILON, Moments, Overlaps, finish, overlap, peak_mean,
};

fn row_of(array: &Array2<f64>, row: usize) -> &[f64] {
    array
        .row(row)
        .to_slice()
        .expect("array row is not contiguous")
}

/// One flat buffer, and the composition half centred once as `f32`. The descent walks pairs in
/// no order it can prefetch, so a vector per row costs a cache miss on the only wide loop.
pub struct PreparedAggregate {
    metric: AggregateMetric,
    aggregation: CoverageAggregation,
    n_samples: usize,
    tnf_width: usize,
    samples: Vec<Moments>,
    presence: Vec<f64>,
    composition: CompositionMetric,
    composition_scale: f64,
    composition_only: bool,
    tnf: Vec<f32>,
    tnf_variance: Vec<f32>,
}

impl PreparedAggregate {
    /// The rows are read out of the two tables rather than a concatenated copy of them, because
    /// the copy is one allocation per contig and the descent throws it away straight after.
    pub fn new(
        coverage_table: &Array2<f64>,
        tnf_table: &Array2<f64>,
        indices: &[usize],
        floors: &[f64],
        settings: DistanceSettings,
    ) -> Self {
        let n_coverage_columns = coverage_table.ncols();
        let n_samples = n_coverage_columns / 2;
        let tnf_width = tnf_table.ncols();
        // At weight zero the arithmetic combination is the composition term alone, so the whole
        // coverage half of the distance is multiplied out and never has to be computed.
        let composition_only = settings.aggregate_weight == Some(0.0)
            && settings.combination == Combination::Arithmetic;

        let held = match composition_only {
            true => 0,
            false => indices.len(),
        };
        let mut samples = Vec::with_capacity(held * n_samples);
        let mut presence = Vec::with_capacity(held);
        let mut tnf = Vec::with_capacity(indices.len() * tnf_width);
        let mut tnf_variance = Vec::with_capacity(indices.len());

        for (index, floor) in indices.iter().zip(floors) {
            if !composition_only {
                let coverage = row_of(coverage_table, *index);
                samples.extend(
                    coverage
                        .chunks_exact(2)
                        .map(|sample| Moments::new(sample[0], (sample[1] + EPSILON).max(*floor))),
                );
                presence.push(settings.presence_fraction * peak_mean(coverage));
            }

            let composition = row_of(tnf_table, *index);
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
            metric: AggregateMetric::new(n_coverage_columns, settings),
            aggregation: settings.aggregation,
            n_samples,
            tnf_width,
            samples,
            presence,
            composition: settings.composition,
            composition_scale: settings.composition_scale,
            composition_only,
            tnf,
            tnf_variance,
        }
    }

    pub fn distance(&self, a: usize, b: usize) -> f64 {
        if self.composition_only {
            let distance = self.composition(a, b);
            return if distance.is_nan() { 1.0 } else { distance };
        }
        let (coverage, scored) = self.coverage(a, b);
        self.metric
            .combine(coverage, scored, self.composition(a, b))
    }

    fn coverage(&self, a: usize, b: usize) -> (f64, usize) {
        let mut overlaps = Overlaps::new(self.aggregation);
        for (x, y) in self.samples_of(a).iter().zip(self.samples_of(b)) {
            if x.mean - EPSILON <= self.presence[a] && y.mean - EPSILON <= self.presence[b] {
                continue;
            }
            overlaps.seen += 1;
            overlaps.push(overlap(*x, *y).clamp(EPSILON, 1.0 - EPSILON));
        }
        finish(&overlaps, false)
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
