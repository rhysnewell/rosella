use ndarray::Array2;

use crate::embedding::features::row_slice;
use crate::embedding::knn::{KnnGraph, Metric};
use crate::learned::device::{Coverage, NOBODY, Rows, nearest};

use super::{
    AggregateMetric, DistanceSettings, EPSILON, MIN_VAR_EPSILON, Moments, coded::Coded,
    coverage_distance, erfc, moments, presence, rho_from,
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
    coded: std::sync::OnceLock<Coded>,
}

impl PreparedAggregate {
    // The rows are read out of the two tables rather than a concatenated copy of them, because
    // the copy is one allocation per contig and the descent throws it away straight after.
    pub fn new(
        coverage_table: &Array2<f64>,
        tnf_table: &Array2<f64>,
        indices: &[usize],
        settings: &DistanceSettings,
    ) -> Self {
        let n_coverage_columns = coverage_table.ncols();
        let n_samples = n_coverage_columns / 2;
        let tnf_width = tnf_table.ncols();
        // At weight zero the combination is the composition term alone, so the whole coverage
        // half of the distance is multiplied out and never has to be computed.
        let composition_only = settings.is_composition_only();

        let held = match composition_only {
            true => 0,
            false => indices.len(),
        };
        let mut prepared = Self {
            metric: AggregateMetric::new(n_coverage_columns, settings.clone()),
            n_samples,
            tnf_width,
            samples: Vec::with_capacity(held * n_samples),
            presence: Vec::with_capacity(held),
            composition_only,
            tnf: Vec::with_capacity(indices.len() * tnf_width),
            tnf_variance: Vec::with_capacity(indices.len()),
            coded: std::sync::OnceLock::new(),
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
                .push(presence(coverage, self.metric.settings().presence_fraction));
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

    // The coded rows bound the composition from below, so a pair already past `bound` there is
    // turned away before its full rows are read.
    fn past(&self, a: usize, b: usize, bound: f64) -> Option<f64> {
        let total = f64::from(self.tnf_variance[a]) + f64::from(self.tnf_variance[b]);
        if !(bound.is_finite() && total > 0.0 && self.tnf_width > 0) {
            return None;
        }
        let coded = self
            .coded
            .get_or_init(|| Coded::of(&self.tnf, self.tnf_width));
        let low = (total - 2.0 * coded.highest_dot(a, b)) / total - 1e-12;
        let low = match self.composition_only {
            true => low,
            false => self.metric.floor(low),
        };
        (low > bound).then_some(low)
    }

    // Queries side by side, so one point is read once and measured against all of them. Each
    // lane still sums in the order `dot` does, so the distances match to the bit.
    pub fn block(&self, queries: &[usize]) -> QueryBlock {
        let mut tnf = vec![0.0f32; self.tnf_width * QUERIES];
        for (lane, query) in queries.iter().take(QUERIES).enumerate() {
            for (at, value) in self.tnf_of(*query).iter().enumerate() {
                tnf[at * QUERIES + lane] = *value;
            }
        }
        QueryBlock {
            queries: queries.iter().take(QUERIES).copied().collect(),
            tnf,
        }
    }

    pub fn within_block(&self, block: &QueryBlock, b: usize, bounds: &[f64], out: &mut [f64]) {
        let mut sums = [-0.0f32; QUERIES];
        for (lanes, y) in block
            .tnf
            .as_chunks::<QUERIES>()
            .0
            .iter()
            .zip(self.tnf_of(b))
        {
            for (sum, x) in sums.iter_mut().zip(lanes) {
                *sum += x * y;
            }
        }
        for (lane, a) in block.queries.iter().enumerate() {
            let composition = self.rho(*a, b, sums[lane]);
            out[lane] = self.settle(*a, b, composition, bounds[lane]);
        }
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

// Composition is a correlation of centred rows, so a device measures every pair in products and
// keeps the exact neighbours. Coverage is added per pair where the composition floor lets it in.
impl PreparedAggregate {
    fn on_device(&self, pairs: usize) -> bool {
        pairs >= DEVICE_PAIRS
    }

    fn rows(&self, end: usize) -> Rows<'_> {
        Rows {
            values: &self.tnf[..end * self.tnf_width],
            variance: &self.tnf_variance[..end],
        }
    }

    // Lanes read neighbouring rows of one sample's field, so the table is laid out field by row.
    fn device_moments(&self) -> Vec<f64> {
        let rows = self.presence.len();
        let mut moments = vec![0.0; rows * self.n_samples * 4];
        for (row, samples) in self.samples.chunks(self.n_samples).enumerate() {
            for (sample, at) in samples.iter().enumerate() {
                let fields = [at.mean, at.variance, at.deviation, at.log_variance];
                for (field, value) in fields.into_iter().enumerate() {
                    moments[(sample * 4 + field) * rows + row] = value;
                }
            }
        }
        moments
    }

    // A negative weight stands for the weight each pair's scored samples give.
    fn device_terms(&self) -> Vec<f64> {
        let settings = self.metric.settings();
        let mut terms = vec![
            settings.aggregate_weight.unwrap_or(-1.0),
            self.metric.heaviest(),
            EPSILON,
            MIN_VAR_EPSILON,
        ];
        terms.extend(erfc::layout());
        terms
    }

    pub fn searched_on_device(&self, base: usize, queries: &[usize], k: usize) -> Option<KnnGraph> {
        let values = queries
            .iter()
            .flat_map(|row| self.tnf_of(*row).iter().copied())
            .collect::<Vec<_>>();
        let variance = queries
            .iter()
            .map(|row| self.tnf_variance[*row])
            .collect::<Vec<_>>();
        let selves = queries
            .iter()
            .map(|row| match *row < base {
                true => *row as u32,
                false => NOBODY,
            })
            .collect::<Vec<_>>();
        let asked = Rows {
            values: &values,
            variance: &variance,
        };
        let rows = queries.iter().map(|row| *row as u32).collect::<Vec<_>>();
        let mixed = (!self.composition_only && self.n_samples > 0)
            .then(|| (self.device_moments(), self.device_terms()));
        let coverage = mixed.as_ref().map(|(moments, terms)| Coverage {
            rows: &rows,
            moments,
            presence: &self.presence,
            samples: self.n_samples,
            terms,
        });
        let (indices, dists) = match nearest(self.rows(base), asked, &selves, coverage, k) {
            Ok(found) => found?,
            Err(error) => {
                log::warn!("GPU: neighbour search failed ({error:#}), searching on the CPU");
                return None;
            }
        };
        let shape = (queries.len(), k);
        Some(KnnGraph {
            indices: Array2::from_shape_vec(shape, indices).ok()?,
            dists: Array2::from_shape_vec(shape, dists).ok()?,
        })
    }

    pub fn device_knn(&self, k: usize) -> Option<KnnGraph> {
        let rows = self.tnf_variance.len();
        if !self.on_device(rows * rows) {
            return None;
        }
        self.searched_on_device(rows, &(0..rows).collect::<Vec<_>>(), k)
    }

    pub fn device_nearest(&self, first: usize, k: usize) -> Option<KnnGraph> {
        let rows = self.tnf_variance.len();
        if !self.on_device(first * (rows - first)) {
            return None;
        }
        self.searched_on_device(first, &(first..rows).collect::<Vec<_>>(), k)
    }

    pub fn device_exact(&self, base: usize, queries: &[usize], k: usize) -> Option<KnnGraph> {
        if !self.on_device(base * queries.len()) {
            return None;
        }
        self.searched_on_device(base, queries, k)
    }

    // A walk reads scattered rows it rarely keeps, so the coded rows turn most away first. Rows
    // already in cache or read in order cost more to filter than to measure.
    fn scattered(
        &self,
        a: usize,
        others: &[u32],
        bound: impl Fn(usize) -> f64,
        out: &mut Vec<f64>,
    ) {
        for chunk in others.chunks(LANES) {
            let mut measured = [0.0; LANES];
            let (mut open, mut lanes, mut count) = ([0; LANES], [0; LANES], 0);
            for (lane, b) in chunk.iter().enumerate() {
                match self.past(a, *b as usize, bound(*b as usize)) {
                    Some(low) => measured[lane] = low,
                    None => {
                        (open[count], lanes[count]) = (*b, lane);
                        count += 1;
                    }
                }
            }
            let compositions = self.compositions(a, &open[..count]);
            for (at, composition) in compositions.iter().take(count).enumerate() {
                let b = open[at] as usize;
                measured[lanes[at]] = self.settle(a, b, *composition, bound(b));
            }
            out.extend_from_slice(&measured[..chunk.len()]);
        }
    }
}

pub struct QueryBlock {
    pub queries: Vec<usize>,
    tnf: Vec<f32>,
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
        self.metric.scattered(a + self.by, others, bound, out);
    }
}

fn dot(a: &[f32], b: &[f32]) -> f32 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

// Each lane sums in the order `dot` does, so the bits match, but the lanes' adds overlap where
// one sum's adds wait on each other.
pub const LANES: usize = 8;

// Under this many pairs the descent finishes before the device has its rows.
const DEVICE_PAIRS: usize = 1 << 30;

// Sixteen sums fill four 128-bit registers, enough to hide the adds' latency without the spills
// thirty-two cost on x86.
pub const QUERIES: usize = 16;

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
