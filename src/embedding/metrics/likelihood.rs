use super::EPSILON;

const TAU: f64 = std::f64::consts::TAU;
const LOG_DENSITY_FLOOR: f64 = -27.631;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum CoverageDistance {
    #[default]
    Overlap,
    Likelihood,
    LikelihoodMean,
    LikelihoodSample,
}

// A log ratio against random pairs rather than a width, so a sample where every contig sits at
// one depth carries little evidence either way.
#[derive(Debug)]
pub struct DepthModel {
    pub form: CoverageDistance,
    pub prior: f64,
    pub samples: Vec<Option<SampleModel>>,
}

#[derive(Debug)]
pub struct SampleModel {
    pub read_length: f64,
    pub spread: f64,
    pub absent: f64,
    pub background: LogDensity,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct LogDepth {
    depth: f64,
    log: f64,
    reads_per_depth: f64,
    variance: f64,
}

impl DepthModel {
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

    pub fn distance(
        &self,
        a: &[LogDepth],
        a_presence: f64,
        b: &[LogDepth],
        b_presence: f64,
    ) -> (f64, usize) {
        let scored = |x: &LogDepth, y: &LogDepth| x.depth > a_presence || y.depth > b_presence;
        let ratios = a
            .iter()
            .zip(b)
            .zip(&self.samples)
            .filter_map(|((x, y), model)| Some((x, y, model.as_ref()?)))
            .filter(|(x, y, _)| scored(x, y))
            .map(|(x, y, model)| model.log_ratio(x, y));
        let (total, count) = match self.form {
            CoverageDistance::LikelihoodSample => {
                let count = ratios.clone().count();
                let share = self.prior / count.max(1) as f64;
                let total = ratios.map(|ratio| logistic(ratio + share)).sum::<f64>();
                (total, count)
            }
            _ => ratios.fold((0.0, 0), |(total, count), ratio| (total + ratio, count + 1)),
        };
        if count == 0 {
            return (EPSILON, 0);
        }
        let distance = match self.form {
            CoverageDistance::LikelihoodSample => total / count as f64,
            CoverageDistance::LikelihoodMean => logistic((total + self.prior) / count as f64),
            _ => logistic(total + self.prior),
        };
        let distance = if distance.is_nan() { 1.0 } else { distance };
        (distance.clamp(EPSILON, 1.0 - EPSILON), count)
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
                let variance = a.variance + b.variance;
                self.background.at(gap) + 0.5 * (gap * gap / variance + (TAU * variance).ln())
            }
            (true, false) => self.absent + a.depth * b.reads_per_depth,
            (false, true) => self.absent + b.depth * a.reads_per_depth,
            (false, false) => 0.0,
        }
    }
}

fn logistic(log_odds: f64) -> f64 {
    1.0 / (1.0 + (-log_odds).exp())
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
