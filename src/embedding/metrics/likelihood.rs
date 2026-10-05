use super::EPSILON;

use std::f64::consts::{PI, TAU};

const LOG_DENSITY_FLOOR: f64 = -27.631;
pub const COMPOSITION_FLOOR: f64 = 1e-6;
pub const TAILS: [f64; 12] = [
    1.0,
    2.0,
    3.0,
    4.0,
    6.0,
    8.0,
    12.0,
    16.0,
    24.0,
    32.0,
    64.0,
    f64::INFINITY,
];

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum CoverageDistance {
    #[default]
    Overlap,
    Likelihood,
}

// Depth and composition each give a log ratio of different over same genome, so they add with no
// fitted share between them.
#[derive(Debug)]
pub struct Evidence {
    pub prior: f64,
    pub samples: Vec<Option<SampleModel>>,
    pub composition: Option<CompositionModel>,
}

#[derive(Debug)]
pub struct SampleModel {
    pub read_length: f64,
    pub spread: f64,
    pub tail: Tail,
    pub absent: f64,
    pub background: LogDensity,
}

// Log composition distance against log(1/a + 1/b) of the two lengths, a line for each class, since
// short contigs scatter wider in both.
#[derive(Debug)]
pub struct CompositionModel {
    pub same: Line,
    pub scale: f64,
    pub tail: Tail,
    pub other: Line,
    pub background: LogDensity,
}

#[derive(Debug, Clone, Copy)]
pub struct Line {
    pub intercept: f64,
    pub slope: f64,
}

impl Line {
    pub fn at(&self, x: f64) -> f64 {
        self.intercept + self.slope * x
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct LogDepth {
    depth: f64,
    log: f64,
    reads_per_depth: f64,
    variance: f64,
}

impl Evidence {
    pub fn depths<'a>(
        &'a self,
        row: &'a [f64],
        length: f64,
    ) -> impl Iterator<Item = LogDepth> + 'a {
        row.chunks_exact(2)
            .zip(&self.samples)
            .map(move |(sample, model)| match model {
                Some(model) => model.depth(sample[0], length),
                None => LogDepth::default(),
            })
    }

    pub fn log_odds(
        &self,
        a: &[LogDepth],
        a_presence: f64,
        b: &[LogDepth],
        b_presence: f64,
    ) -> (f64, usize) {
        a.iter()
            .zip(b)
            .zip(&self.samples)
            .filter_map(|((x, y), model)| Some((x, y, model.as_ref()?)))
            .filter(|(x, y, _)| x.depth > a_presence || y.depth > b_presence)
            .fold((0.0, 0), |(total, count), (x, y, model)| {
                (total + model.log_ratio(x, y), count + 1)
            })
    }

    pub fn composition_ratio(&self, distance: f64, a_length: f64, b_length: f64) -> f64 {
        self.composition
            .as_ref()
            .map_or(0.0, |model| model.log_ratio(distance, a_length, b_length))
    }

    pub fn posterior(&self, log_odds: f64) -> f64 {
        let distance = 1.0 / (1.0 + (-(log_odds + self.prior)).exp());
        if distance.is_nan() { 1.0 } else { distance }
    }

    pub fn coverage(&self, log_odds: f64, scored: usize) -> (f64, usize) {
        match scored {
            0 => (EPSILON, 0),
            _ => (
                self.posterior(log_odds).clamp(EPSILON, 1.0 - EPSILON),
                scored,
            ),
        }
    }
}

impl SampleModel {
    fn depth(&self, depth: f64, length: f64) -> LogDepth {
        let reads_per_depth = length / self.read_length;
        match depth > 0.0 {
            true => LogDepth {
                depth,
                log: depth.ln(),
                reads_per_depth,
                variance: 1.0 / (depth * reads_per_depth) + self.spread,
            },
            false => LogDepth {
                reads_per_depth,
                ..LogDepth::default()
            },
        }
    }

    // A zero has no log, so it is scored as a read count against the share of contigs with none.
    fn log_ratio(&self, a: &LogDepth, b: &LogDepth) -> f64 {
        match (a.depth > 0.0, b.depth > 0.0) {
            (true, true) => {
                let gap = a.log - b.log;
                self.background.at(gap) - self.tail.log_density(gap, a.variance + b.variance)
            }
            (true, false) => self.absent + a.depth * b.reads_per_depth,
            (false, true) => self.absent + b.depth * a.reads_per_depth,
            (false, false) => 0.0,
        }
    }
}

impl CompositionModel {
    pub fn log_ratio(&self, distance: f64, a_length: f64, b_length: f64) -> f64 {
        let x = distance.max(COMPOSITION_FLOOR).ln();
        let u = length_term(a_length, b_length);
        self.background.at(x - self.other.at(u))
            - self.tail.log_density(x - self.same.at(u), self.scale)
    }
}

pub fn length_term(a_length: f64, b_length: f64) -> f64 {
    (1.0 / a_length + 1.0 / b_length).ln()
}

// Genome mates stray far more often than a normal allows, and a normal makes one stray sample
// outvote every other, so the same-genome model is a Student t whose tail the data picks.
#[derive(Debug, Clone, Copy)]
pub struct Tail {
    pub degrees: f64,
    constant: f64,
}

impl Tail {
    pub fn new(degrees: f64) -> Self {
        let constant = match degrees.is_finite() {
            true => log_gamma_ratio(degrees) - 0.5 * (degrees * PI).ln(),
            false => -0.5 * TAU.ln(),
        };
        Self { degrees, constant }
    }

    pub fn log_density(&self, x: f64, scale: f64) -> f64 {
        let spread = x * x / scale;
        match self.degrees.is_finite() {
            true => {
                self.constant
                    - 0.5 * scale.ln()
                    - 0.5 * (self.degrees + 1.0) * (spread / self.degrees).ln_1p()
            }
            false => self.constant - 0.5 * (scale.ln() + spread),
        }
    }
}

// ln Gamma((v + 1) / 2) - ln Gamma(v / 2) for whole v, stepped up from v = 1 or v = 2.
fn log_gamma_ratio(degrees: f64) -> f64 {
    let whole = degrees as u32;
    let (mut at, mut ratio) = match whole % 2 {
        1 => (1, -0.5 * PI.ln()),
        _ => (2, 0.5 * PI.ln() - 2f64.ln()),
    };
    while at < whole {
        ratio += ((at + 1) as f64 / at as f64).ln();
        at += 2;
    }
    ratio
}

#[derive(Debug)]
pub struct LogDensity {
    pub start: f64,
    pub step: f64,
    pub values: Vec<f64>,
}

impl LogDensity {
    pub fn of(start: f64, step: f64, density: &[f64]) -> Self {
        Self {
            start,
            step,
            values: density
                .iter()
                .map(|value| value.ln().max(LOG_DENSITY_FLOOR))
                .collect(),
        }
    }

    pub fn at(&self, x: f64) -> f64 {
        let at = (x - self.start) / self.step;
        if !(0.0..(self.values.len() - 1) as f64).contains(&at) {
            return LOG_DENSITY_FLOOR;
        }
        let left = at.floor() as usize;
        let right = at - left as f64;
        self.values[left] * (1.0 - right) + self.values[left + 1] * right
    }
}
