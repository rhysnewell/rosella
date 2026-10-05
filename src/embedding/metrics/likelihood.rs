use super::EPSILON;

const TAU: f64 = std::f64::consts::TAU;
const LOG_DENSITY_FLOOR: f64 = -27.631;
pub const COMPOSITION_FLOOR: f64 = 1e-6;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum CoverageDistance {
    #[default]
    Overlap,
    Likelihood,
}

// Depth and composition each give a log ratio of different over same genome, so they add with no
// fitted share between them.
#[derive(Debug, Clone)]
pub struct Evidence {
    pub prior: f64,
    pub samples: Vec<Option<SampleModel>>,
    pub composition: Option<CompositionModel>,
}

#[derive(Debug, Clone)]
pub struct SampleModel {
    pub read_length: f64,
    pub spread: f64,
    pub absent: f64,
    pub background: LogDensity,
}

// Log composition distance against log(1/a + 1/b) of the two lengths, a line for each class, since
// short contigs scatter wider in both.
#[derive(Debug, Clone)]
pub struct CompositionModel {
    pub same: Line,
    pub spread: f64,
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
                self.background.at(gap) - log_normal(gap, a.variance + b.variance)
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
            - log_normal(x - self.same.at(u), self.spread * self.spread)
    }
}

pub fn length_term(a_length: f64, b_length: f64) -> f64 {
    (1.0 / a_length + 1.0 / b_length).ln()
}

pub fn log_normal(x: f64, variance: f64) -> f64 {
    -0.5 * (x * x / variance + (TAU * variance).ln())
}

#[derive(Debug, Clone)]
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
