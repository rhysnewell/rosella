use crate::learned::labelled::net::Params;

const SECOND_MOMENT: f64 = 0.999;
const EPSILON: f64 = 1e-8;

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
        }
    }
}
