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
}

fn cosine(from: f64, to: f64, share: f64) -> f64 {
    to + (from - to) / 2.0 * (1.0 + (PI * share).cos())
}

impl Schedule {
    pub fn new(steps: usize, peak: f64) -> Self {
        Self { steps, peak }
    }

    pub fn at(&self, step: usize) -> (f64, f64) {
        let start = self.peak / START_DIVISOR;
        let end = start / END_DIVISOR;
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

// Decoupled weight decay (Loshchilov and Hutter 2019).
pub struct AdamW {
    first: Params,
    second: Params,
    step: i32,
    decay: f64,
}

impl AdamW {
    pub fn new(params: &Params, decay: f64) -> Self {
        Self {
            first: params.zeros_like(),
            second: params.zeros_like(),
            step: 0,
            decay,
        }
    }

    pub fn update(&mut self, params: &mut Params, grads: &Params, rate: f64, momentum: f64) {
        self.step += 1;
        let first_correction = 1.0 - momentum.powi(self.step);
        let second_correction = (1.0 - SECOND_MOMENT.powi(self.step)).sqrt();
        let shrink = (1.0 - rate * self.decay) as f32;
        let step = (rate / first_correction) as f32;
        let (momentum, epsilon) = (momentum as f32, EPSILON as f32);
        let (keep, correction) = (SECOND_MOMENT as f32, second_correction as f32);
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
