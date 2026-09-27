use rayon::prelude::*;

const GRID: usize = 64;
const DENSITY_FLOOR: f64 = 1e-12;
const EM_ROUNDS: usize = 500;
const NEWTON_STEPS: usize = 5;
const SETTLED: f64 = 1e-7;

#[derive(Clone, Copy, Debug)]
pub struct Sample {
    pub share: f64,
    pub length: usize,
}

pub struct Fit {
    pub chances: Vec<f64>,
    pub prior: Vec<(usize, f64)>,
}

// Pieces searched with and without their own bin show a contig with and without a home at each
// length, so the real contigs' mix of the two is fitted on this assembly and nothing is carried in.
pub fn chances(homed: &[(Sample, bool)], homeless: &[Sample], real: &[Sample]) -> Option<Fit> {
    if homed.is_empty() || homeless.is_empty() || real.is_empty() {
        return None;
    }
    let log_length = |sample: &Sample| (sample.length.max(1) as f64).ln();
    let (low, high) = homed
        .iter()
        .map(|(sample, _)| sample)
        .chain(homeless)
        .chain(real)
        .map(log_length)
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(low, high), y| {
            (low.min(y), high.max(y))
        });
    let grid = Grid::new(low, high);
    let points = |samples: &mut dyn Iterator<Item = &Sample>| {
        samples
            .map(|sample| (sample.share, log_length(sample)))
            .collect::<Vec<_>>()
    };
    let homed_points = points(&mut homed.iter().map(|(sample, _)| sample));
    let found = homed
        .iter()
        .map(|(_, right)| f64::from(u8::from(*right)))
        .collect::<Vec<_>>();
    let with_home = grid.density(&homed_points, None)?;
    let right_mass = grid.density(&homed_points, Some(&found))?;
    let without_home = grid.density(&points(&mut homeless.iter()), None)?;

    let real_points = points(&mut real.iter());
    let at = real_points
        .iter()
        .map(|(x, y)| {
            let home = grid.at(&with_home, *x, *y);
            let right = match home > 0.0 {
                true => (grid.at(&right_mass, *x, *y) / home).clamp(0.0, 1.0),
                false => 0.0,
            };
            (
                home.max(DENSITY_FLOOR),
                grid.at(&without_home, *x, *y).max(DENSITY_FLOOR),
                right,
            )
        })
        .collect::<Vec<_>>();
    let centre = real_points.iter().map(|(_, y)| y).sum::<f64>() / real_points.len() as f64;
    let offsets = real_points
        .iter()
        .map(|(_, y)| y - centre)
        .collect::<Vec<_>>();
    let (intercept, slope) = prior_by_length(&at, &offsets);
    let prior = |offset: f64| logistic(intercept + slope * offset);
    let chances = at
        .iter()
        .zip(&offsets)
        .map(|((home, away, right), offset)| {
            let p = prior(*offset);
            p * home / (p * home + (1.0 - p) * away) * right
        })
        .collect();
    let mut lengths = real.iter().map(|sample| sample.length).collect::<Vec<_>>();
    lengths.sort_unstable();
    let prior = [0.1, 0.5, 0.9]
        .iter()
        .map(|quantile| {
            let length = lengths[((lengths.len() - 1) as f64 * quantile) as usize];
            (length, prior((length.max(1) as f64).ln() - centre))
        })
        .collect();
    Some(Fit { chances, prior })
}

fn prior_by_length(at: &[(f64, f64, f64)], offsets: &[f64]) -> (f64, f64) {
    let (mut intercept, mut slope) = (0.0, 0.0);
    let mut held = vec![0.5; at.len()];
    for _ in 0..EM_ROUNDS {
        let weights = at
            .iter()
            .zip(offsets)
            .map(|((home, away, _), offset)| {
                let p = logistic(intercept + slope * offset);
                p * home / (p * home + (1.0 - p) * away)
            })
            .collect::<Vec<_>>();
        for _ in 0..NEWTON_STEPS {
            let (mut g0, mut g1, mut h00, mut h01, mut h11) = (0.0, 0.0, 0.0, 0.0, 0.0);
            for (weight, offset) in weights.iter().zip(offsets) {
                let p = logistic(intercept + slope * offset);
                let residual = weight - p;
                let curvature = p * (1.0 - p);
                g0 += residual;
                g1 += residual * offset;
                h00 += curvature;
                h01 += curvature * offset;
                h11 += curvature * offset * offset;
            }
            let determinant = h00 * h11 - h01 * h01;
            if determinant.abs() < f64::EPSILON {
                break;
            }
            intercept += (h11 * g0 - h01 * g1) / determinant;
            slope += (h00 * g1 - h01 * g0) / determinant;
        }
        let moved = weights
            .iter()
            .zip(&held)
            .map(|(now, before)| (now - before).abs())
            .fold(0.0, f64::max);
        held = weights;
        if moved < SETTLED {
            break;
        }
    }
    (intercept, slope)
}

fn logistic(value: f64) -> f64 {
    1.0 / (1.0 + (-value.clamp(-50.0, 50.0)).exp())
}

struct Grid {
    low: f64,
    step_y: f64,
}

impl Grid {
    fn new(low: f64, high: f64) -> Self {
        Self {
            low,
            step_y: ((high - low) / (GRID - 1) as f64).max(f64::EPSILON),
        }
    }

    fn node(&self, at: usize) -> (f64, f64) {
        let (i, j) = (at / GRID, at % GRID);
        (
            i as f64 / (GRID - 1) as f64,
            self.low + j as f64 * self.step_y,
        )
    }

    // Never narrower than one grid step, so a column of identical lengths still spreads.
    fn density(&self, points: &[(f64, f64)], weights: Option<&[f64]>) -> Option<Vec<f64>> {
        let n = points.len();
        if n < 2 {
            return None;
        }
        let spread = |values: &mut dyn Iterator<Item = f64>| {
            let values = values.collect::<Vec<_>>();
            let mean = values.iter().sum::<f64>() / n as f64;
            (values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1) as f64).sqrt()
        };
        let scale = (n as f64).powf(-1.0 / 6.0);
        let hx = (spread(&mut points.iter().map(|p| p.0)) * scale).max(1.0 / (GRID - 1) as f64);
        let hy = (spread(&mut points.iter().map(|p| p.1)) * scale).max(self.step_y);
        let norm = 1.0 / (n as f64 * 2.0 * std::f64::consts::PI * hx * hy);
        Some(
            (0..GRID * GRID)
                .into_par_iter()
                .map(|at| {
                    let (x, y) = self.node(at);
                    points
                        .iter()
                        .enumerate()
                        .map(|(k, (px, py))| {
                            let weight = weights.map_or(1.0, |weights| weights[k]);
                            weight
                                * (-0.5 * (((x - px) / hx).powi(2) + ((y - py) / hy).powi(2))).exp()
                        })
                        .sum::<f64>()
                        * norm
                })
                .collect(),
        )
    }

    fn at(&self, values: &[f64], x: f64, y: f64) -> f64 {
        let fx = (x.clamp(0.0, 1.0) * (GRID - 1) as f64).min((GRID - 1) as f64 - 1e-9);
        let fy = ((y - self.low) / self.step_y).clamp(0.0, (GRID - 1) as f64 - 1e-9);
        let (i, j) = (fx.floor() as usize, fy.floor() as usize);
        let (dx, dy) = (fx - i as f64, fy - j as f64);
        let value = |i: usize, j: usize| values[i * GRID + j];
        let next_j = (j + 1).min(GRID - 1);
        let next_i = (i + 1).min(GRID - 1);
        value(i, j) * (1.0 - dx) * (1.0 - dy)
            + value(next_i, j) * dx * (1.0 - dy)
            + value(i, next_j) * (1.0 - dx) * dy
            + value(next_i, next_j) * dx * dy
    }
}
