use std::f64::consts::PI;

use rayon::prelude::*;

use crate::learned::net::Params;

const SECOND_MOMENT: f64 = 0.999;
const EPSILON: f64 = 1e-8;
const WARM: f64 = 0.3;
const START_DIVISOR: f64 = 25.0;
const END_DIVISOR: f64 = 1e4;
const MOMENTA: (f64, f64) = (0.85, 0.95);
const CHUNK: usize = 1 << 14;

// One cycle (Smith and Topin 2019): the rate climbs over the first 0.3 of the steps and anneals
// over the rest, with momentum moving the other way, as PyTorch's OneCycleLR runs it.
pub struct Schedule {
    steps: usize,
    peak: f64,
    tail: Option<Tail>,
}

struct Tail {
    from: usize,
    steps: usize,
    rate: f64,
    momentum: f64,
}

fn cosine(from: f64, to: f64, share: f64) -> f64 {
    to + (from - to) / 2.0 * (1.0 + (PI * share).cos())
}

impl Schedule {
    pub fn new(steps: usize, peak: f64) -> Self {
        Self {
            steps,
            peak,
            tail: None,
        }
    }

    pub fn annealing(&self, step: usize) -> bool {
        step as f64 > WARM * self.steps as f64 - 1.0
    }

    // A run that stops early still anneals, from wherever the rate stands down to the end rate.
    pub fn anneal_from(&mut self, step: usize, steps: usize) {
        let (rate, momentum) = self.at(step);
        self.tail = Some(Tail {
            from: step,
            steps,
            rate,
            momentum,
        });
    }

    pub fn at(&self, step: usize) -> (f64, f64) {
        let start = self.peak / START_DIVISOR;
        let end = start / END_DIVISOR;
        if let Some(tail) = self.tail.as_ref().filter(|tail| step >= tail.from) {
            let share =
                ((step - tail.from) as f64 / tail.steps.saturating_sub(1).max(1) as f64).min(1.0);
            return (
                cosine(tail.rate, end, share),
                cosine(tail.momentum, MOMENTA.1, share),
            );
        }
        let turn = WARM * self.steps as f64 - 1.0;
        let last = self.steps as f64 - 1.0;
        let step = step as f64;
        if step <= turn {
            let share = if turn > 0.0 { step / turn } else { 1.0 };
            (
                cosine(start, self.peak, share),
                cosine(MOMENTA.1, MOMENTA.0, share),
            )
        } else {
            let share = (step - turn) / (last - turn).max(f64::EPSILON);
            (
                cosine(self.peak, end, share),
                cosine(MOMENTA.0, MOMENTA.1, share),
            )
        }
    }
}

pub fn epoch_message(epoch: usize, epochs: usize, loss: Option<f64>) -> String {
    match loss {
        Some(loss) => format!("epoch {}/{epochs}, loss {loss:.4}", epoch + 1),
        None => format!("epoch {}/{epochs}", epoch + 1),
    }
}

// The noise is the scatter of the last annealed epoch means about their line. Step losses are no
// guide to it. Every epoch repeats their spread over batch sizes and contig mixes.
pub struct Plateau {
    sum: f64,
    steps: usize,
    means: Vec<f64>,
    fall: Option<f64>,
    patience: usize,
}

impl Plateau {
    pub fn new(patience: usize) -> Self {
        Self {
            sum: 0.0,
            steps: 0,
            means: Vec::new(),
            fall: None,
            patience,
        }
    }

    // A drawn run stops with the rate still near its peak, so it anneals over the patience window.
    pub fn tail(&self, drawn: bool) -> usize {
        match drawn {
            true => self.patience.max(1),
            false => 1,
        }
    }

    pub fn shown(&self) -> String {
        self.fall
            .map_or_else(String::new, |fall| format!(" z:{fall:.1}"))
    }

    pub fn record(&mut self, loss: f64) {
        self.sum += loss;
        self.steps += 1;
    }

    pub fn close(&mut self, annealed: bool) -> (f64, bool) {
        let mean = self.sum / self.steps.max(1) as f64;
        (self.sum, self.steps) = (0.0, 0);
        if self.patience == 0 || !annealed {
            return (mean, false);
        }
        self.means.push(mean);
        if self.means.len() <= self.patience {
            return (mean, false);
        }
        let window = &self.means[self.means.len() - self.patience - 1..];
        let (slope, error) = fit(window);
        let fall = -slope / error;
        self.fall = Some(fall);
        (mean, fall.is_nan() || fall < 2.0)
    }
}

fn fit(window: &[f64]) -> (f64, f64) {
    let count = window.len() as f64;
    let centre = (count - 1.0) / 2.0;
    let level = window.iter().sum::<f64>() / count;
    let spread = (0..window.len())
        .map(|at| (at as f64 - centre).powi(2))
        .sum::<f64>();
    let slope = window
        .iter()
        .enumerate()
        .map(|(at, mean)| (at as f64 - centre) * (mean - level))
        .sum::<f64>()
        / spread;
    let scatter = window
        .iter()
        .enumerate()
        .map(|(at, mean)| (mean - level - slope * (at as f64 - centre)).powi(2))
        .sum::<f64>();
    (slope, (scatter / (count - 2.0).max(1.0) / spread).sqrt())
}

pub(crate) struct Factors {
    pub(crate) shrink: f32,
    pub(crate) step: f32,
    pub(crate) momentum: f32,
    pub(crate) keep: f32,
    pub(crate) correction: f32,
    pub(crate) epsilon: f32,
}

pub(crate) struct Clock {
    step: i32,
    decay: f64,
}

impl Clock {
    pub(crate) fn new(decay: f64) -> Self {
        Self { step: 0, decay }
    }

    pub(crate) fn tick(&mut self, rate: f64, momentum: f64) -> Factors {
        self.step += 1;
        let first_correction = 1.0 - momentum.powi(self.step);
        let second_correction = (1.0 - SECOND_MOMENT.powi(self.step)).sqrt();
        Factors {
            shrink: (1.0 - rate * self.decay) as f32,
            step: (rate / first_correction) as f32,
            momentum: momentum as f32,
            keep: SECOND_MOMENT as f32,
            correction: second_correction as f32,
            epsilon: EPSILON as f32,
        }
    }
}

// Decoupled weight decay (Loshchilov and Hutter 2019).
pub struct AdamW {
    first: Params,
    second: Params,
    clock: Clock,
}

impl AdamW {
    pub fn new(params: &Params, decay: f64) -> Self {
        Self {
            first: params.zeros_like(),
            second: params.zeros_like(),
            clock: Clock::new(decay),
        }
    }

    pub fn update(&mut self, params: &mut Params, grads: &Params, rate: f64, momentum: f64) {
        let Factors {
            shrink,
            step,
            momentum,
            keep,
            correction,
            epsilon,
        } = self.clock.tick(rate, momentum);
        for (((params, grads), first), second) in params
            .slices_mut()
            .into_iter()
            .zip(grads.slices())
            .zip(self.first.slices_mut())
            .zip(self.second.slices_mut())
        {
            params
                .par_chunks_mut(CHUNK)
                .zip(grads.par_chunks(CHUNK))
                .zip(first.par_chunks_mut(CHUNK))
                .zip(second.par_chunks_mut(CHUNK))
                .for_each(|(((params, grads), first), second)| {
                    for (((param, grad), first), second) in params
                        .iter_mut()
                        .zip(grads)
                        .zip(first.iter_mut())
                        .zip(second.iter_mut())
                    {
                        *param *= shrink;
                        *first = momentum * *first + (1.0 - momentum) * grad;
                        *second = keep * *second + (1.0 - keep) * grad * grad;
                        *param -= step * *first / (second.sqrt() / correction + epsilon);
                    }
                });
        }
    }
}
