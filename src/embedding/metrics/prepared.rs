use ndarray::Array2;

use crate::embedding::features::row_slice;
use crate::embedding::knn::Metric;

use super::{
    AggregateMetric, DistanceSettings, Moments, coverage_distance, moments, presence, rho_from,
};

// One flat buffer, and the composition half centred once as `f32`. The descent walks pairs in
// no order it can prefetch, so a vector per row costs a cache miss on the only wide loop.
pub struct PreparedAggregate {
    metric: AggregateMetric,
    n_samples: usize,
    tnf_width: usize,
    samples: Vec<Moments>,
    presence: Vec<f64>,
    composition_only: bool,
    tnf: Vec<f32>,
    tnf_variance: Vec<f32>,
}

impl PreparedAggregate {
    // The rows are read out of the two tables rather than a concatenated copy of them, because
    // the copy is one allocation per contig and the descent throws it away straight after.
    pub fn new(
        coverage_table: &Array2<f64>,
        tnf_table: &Array2<f64>,
        indices: &[usize],
        settings: DistanceSettings,
    ) -> Self {
        let n_coverage_columns = coverage_table.ncols();
        let n_samples = n_coverage_columns / 2;
        let tnf_width = tnf_table.ncols();
        // At weight zero the combination is the composition term alone, so the whole coverage
        // half of the distance is multiplied out and never has to be computed.
        let composition_only = settings.aggregate_weight == Some(0.0);

        let held = match composition_only {
            true => 0,
            false => indices.len(),
        };
        let mut prepared = Self {
            metric: AggregateMetric::new(n_coverage_columns, settings),
            n_samples,
            tnf_width,
            samples: Vec::with_capacity(held * n_samples),
            presence: Vec::with_capacity(held),
            composition_only,
            tnf: Vec::with_capacity(indices.len() * tnf_width),
            tnf_variance: Vec::with_capacity(indices.len()),
        };
        for index in indices {
            prepared.push(
                row_slice(coverage_table, *index),
                row_slice(tnf_table, *index),
            );
        }
        prepared
    }

    fn push(&mut self, coverage: &[f64], composition: &[f64]) {
        if !self.composition_only {
            self.samples.extend(moments(coverage));
            self.presence
                .push(presence(coverage, self.metric.presence_fraction()));
        }
        let mean = match composition.is_empty() {
            true => 0.0,
            false => composition.iter().sum::<f64>() / composition.len() as f64,
        };
        let start = self.tnf.len();
        self.tnf
            .extend(composition.iter().map(|value| (value - mean) as f32));
        self.tnf_variance
            .push(dot(&self.tnf[start..], &self.tnf[start..]));
    }

    pub fn distance(&self, a: usize, b: usize) -> f64 {
        self.settle(a, b, self.composition(a, b), f64::INFINITY)
    }

    // Composition costs a dot product and coverage an erfc per sample, so the cheap half
    // decides first whether the dear one can matter.
    fn settle(&self, a: usize, b: usize, composition: f64, bound: f64) -> f64 {
        if self.composition_only {
            return if composition.is_nan() {
                1.0
            } else {
                composition
            };
        }
        let floor = self.metric.floor(composition);
        if floor > bound {
            return floor;
        }
        let (coverage, scored) = self.coverage(a, b);
        self.metric.combine(coverage, scored, composition)
    }

    pub fn shifted(&self, by: usize) -> Shifted<'_> {
        Shifted { metric: self, by }
    }

    fn coverage(&self, a: usize, b: usize) -> (f64, usize) {
        coverage_distance(
            self.samples_of(a),
            self.presence[a],
            self.samples_of(b),
            self.presence[b],
        )
    }

    fn composition(&self, a: usize, b: usize) -> f64 {
        self.rho(a, b, dot(self.tnf_of(a), self.tnf_of(b)))
    }

    fn rho(&self, a: usize, b: usize, dot: f32) -> f64 {
        rho_from(
            dot as f64,
            self.tnf_variance[a] as f64,
            self.tnf_variance[b] as f64,
        )
    }

    fn compositions(&self, a: usize, others: &[u32]) -> [f64; LANES] {
        let mine = self.tnf_of(a);
        let dots = match <&[u32; LANES]>::try_from(others) {
            Ok(full) => dots(mine, full.map(|b| self.tnf_of(b as usize))),
            Err(_) => std::array::from_fn(|lane| {
                others
                    .get(lane)
                    .map_or(0.0, |b| dot(mine, self.tnf_of(*b as usize)))
            }),
        };
        let mut compositions = [0.0; LANES];
        for ((composition, b), dot) in compositions.iter_mut().zip(others).zip(dots) {
            let b = *b as usize;
            *composition = self.rho(a, b, dot);
        }
        compositions
    }

    fn samples_of(&self, row: usize) -> &[Moments] {
        &self.samples[row * self.n_samples..(row + 1) * self.n_samples]
    }

    fn tnf_of(&self, row: usize) -> &[f32] {
        &self.tnf[row * self.tnf_width..(row + 1) * self.tnf_width]
    }
}

impl Metric for &PreparedAggregate {
    fn distance(&self, a: usize, b: usize) -> f64 {
        PreparedAggregate::distance(self, a, b)
    }

    fn within(&self, a: usize, b: usize, bound: f64) -> f64 {
        self.settle(a, b, self.composition(a, b), bound)
    }

    fn within_many(
        &self,
        a: usize,
        others: &[u32],
        bound: impl Fn(usize) -> f64,
        out: &mut Vec<f64>,
    ) {
        for chunk in others.chunks(LANES) {
            let compositions = self.compositions(a, chunk);
            for (b, composition) in chunk.iter().zip(compositions) {
                let b = *b as usize;
                out.push(self.settle(a, b, composition, bound(b)));
            }
        }
    }
}

pub struct Shifted<'a> {
    metric: &'a PreparedAggregate,
    by: usize,
}

impl Metric for Shifted<'_> {
    fn distance(&self, a: usize, b: usize) -> f64 {
        self.metric.distance(a + self.by, b)
    }

    fn within(&self, a: usize, b: usize, bound: f64) -> f64 {
        self.metric.within(a + self.by, b, bound)
    }

    fn within_many(
        &self,
        a: usize,
        others: &[u32],
        bound: impl Fn(usize) -> f64,
        out: &mut Vec<f64>,
    ) {
        self.metric.within_many(a + self.by, others, bound, out);
    }
}

fn dot(a: &[f32], b: &[f32]) -> f32 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

// Each lane sums in the order `dot` does, so the bits match, but the lanes' adds overlap where
// one sum's adds wait on each other.
const LANES: usize = 8;

fn dots(a: &[f32], rows: [&[f32]; LANES]) -> [f32; LANES] {
    let rows = rows.map(|row| &row[..a.len()]);
    let mut sums = [-0.0f32; LANES];
    for (at, x) in a.iter().enumerate() {
        for (sum, row) in sums.iter_mut().zip(&rows) {
            *sum += x * row[at];
        }
    }
    sums
}
