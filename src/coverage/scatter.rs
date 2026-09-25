use ndarray::Array2;
use rayon::prelude::*;

mod pairs;

pub use pairs::fit_neighbours;

const LEAST_MEMBERS: usize = 5;
const LEAST_POINTS: usize = 200;
const MOST_POINTS: usize = 20_000;
const STEPS: usize = 41;
const SAMPLING_EXPONENTS: (f64, f64) = (0.0, 5.0);
const BIAS_EXPONENTS: (f64, f64) = (-6.0, 0.0);
// Heavy tails so the passengers a near complete bin still carries do not widen the fit.
const FREEDOM: f64 = 4.0;

#[derive(Debug, Clone, PartialEq)]
pub enum Source {
    Bins,
    Composition,
    Neighbours,
    Given(Vec<Scatter>),
}

impl Source {
    pub fn parse(choice: &str, given: Option<&str>) -> anyhow::Result<Option<Self>> {
        if let Some(given) = given {
            return Ok(Some(Self::Given(
                given
                    .split(',')
                    .map(Scatter::parse)
                    .collect::<anyhow::Result<_>>()?,
            )));
        }
        Ok(match choice {
            "bins" => Some(Self::Bins),
            "composition" => Some(Self::Composition),
            "neighbours" => Some(Self::Neighbours),
            _ => None,
        })
    }
}

// Read sampling shrinks with contig length and bias along the genome grows with depth, so a
// contig's depth strays from its genome's by both.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Scatter {
    pub sampling: f64,
    pub bias: f64,
}

impl Scatter {
    fn parse(text: &str) -> anyhow::Result<Self> {
        let (sampling, bias) = text
            .split_once(':')
            .ok_or_else(|| anyhow::anyhow!("{text} is not sampling:bias"))?;
        Ok(Self {
            sampling: sampling.trim().parse()?,
            bias: bias.trim().parse()?,
        })
    }

    pub fn variance(&self, depth: f64, length: usize) -> f64 {
        let length = length.max(1) as f64;
        let depth = depth.max(0.0) + self.sampling / length;
        self.sampling * depth / length + self.bias * depth * depth
    }
}

pub fn fit(table: &Array2<f64>, lengths: &[usize], bins: &[Vec<usize>]) -> Vec<Option<Scatter>> {
    (0..table.ncols() / 2)
        .map(|sample| fit_sample(table, lengths, bins, sample))
        .collect()
}

// A sample with no fit keeps today's floor, so a zero variance column never turns sharp.
pub fn apply(table: &mut Array2<f64>, lengths: &[usize], models: &[Option<Scatter>], floor: f64) {
    for (mut row, length) in table.rows_mut().into_iter().zip(lengths) {
        for (sample, model) in models.iter().enumerate() {
            row[2 * sample + 1] = match model {
                Some(model) => model.variance(row[2 * sample], *length),
                None => row[2 * sample + 1].max(floor),
            };
        }
    }
}

struct Point {
    centre: f64,
    length: usize,
    residual: f64,
}

fn fit_sample(
    table: &Array2<f64>,
    lengths: &[usize],
    bins: &[Vec<usize>],
    sample: usize,
) -> Option<Scatter> {
    let depth = |contig: usize| table[[contig, 2 * sample]];
    let mut points = Vec::new();
    for bin in bins.iter().filter(|bin| bin.len() >= LEAST_MEMBERS) {
        let centre = median(bin.iter().map(|contig| depth(*contig)).collect());
        if centre.is_nan() || centre <= 0.0 {
            continue;
        }
        points.extend(bin.iter().map(|contig| Point {
            centre,
            length: lengths[*contig],
            residual: depth(*contig) - centre,
        }));
    }
    if points.len() < LEAST_POINTS {
        return None;
    }
    let stride = points.len().div_ceil(MOST_POINTS);
    let points = points.into_iter().step_by(stride).collect::<Vec<_>>();
    candidates(STEPS)
        .into_par_iter()
        .map(|model| (loss(&model, &points), model))
        .min_by(|a, b| a.0.total_cmp(&b.0))
        .map(|(_, model)| model)
}

fn candidates(steps: usize) -> Vec<Scatter> {
    grid(SAMPLING_EXPONENTS, steps)
        .flat_map(|sampling| {
            grid(BIAS_EXPONENTS, steps).map(move |bias| Scatter { sampling, bias })
        })
        .filter(|model| model.sampling > 0.0 || model.bias > 0.0)
        .collect()
}

fn grid(exponents: (f64, f64), steps: usize) -> impl Iterator<Item = f64> {
    std::iter::once(0.0).chain((0..steps).map(move |step| {
        10f64.powf(exponents.0 + (exponents.1 - exponents.0) * step as f64 / (steps - 1) as f64)
    }))
}

fn loss(model: &Scatter, points: &[Point]) -> f64 {
    points
        .iter()
        .map(|point| {
            let variance = model.sampling * point.centre / point.length.max(1) as f64
                + model.bias * point.centre * point.centre;
            0.5 * variance.ln()
                + 0.5
                    * (FREEDOM + 1.0)
                    * (point.residual * point.residual / (FREEDOM * variance)).ln_1p()
        })
        .sum()
}

fn median(mut values: Vec<f64>) -> f64 {
    let middle = values.len() / 2;
    let (_, value, _) = values.select_nth_unstable_by(middle, f64::total_cmp);
    *value
}
