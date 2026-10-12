use crate::learned::gemm::{fresh, later, later_turned, product};
use crate::learned::norm::{self, Normed, Running};
use ndarray::{Array1, Array2, ArrayView2, s};
use rand::{Rng, rngs::StdRng};

pub const WIDTH: usize = 128;
// Each layer halves the one before. Against a square 2048 second layer this read higher on every
// tier of mouse gut and trained faster.
const LAYERS: [usize; 4] = [2048, 1024, 512, 256];

#[derive(Clone)]
struct Layer {
    weights: Array2<f32>,
    scale: Array1<f32>,
    shift: Array1<f32>,
}

#[derive(Clone)]
pub struct Params {
    layers: Vec<Layer>,
    last: Array2<f32>,
    last_bias: Array1<f32>,
}

impl Params {
    // The uniform bound PyTorch gives a linear layer by default.
    fn new(input: usize, rng: &mut StdRng) -> Self {
        let mut uniform = |rows: usize, columns: usize| {
            let bound = 1.0 / (columns as f32).sqrt();
            Array2::from_shape_simple_fn((rows, columns), || rng.random_range(-bound..bound))
        };
        let mut from = input;
        let mut layers = Vec::with_capacity(LAYERS.len());
        for width in LAYERS {
            layers.push(Layer {
                weights: uniform(width, from),
                scale: Array1::ones(width),
                shift: Array1::zeros(width),
            });
            from = width;
        }
        let last = uniform(WIDTH, from);
        let last_bias = uniform(1, from).slice(s![0, ..WIDTH]).to_owned();
        Self {
            layers,
            last,
            last_bias,
        }
    }

    pub fn zeros_like(&self) -> Self {
        Self {
            layers: self
                .layers
                .iter()
                .map(|layer| Layer {
                    weights: Array2::zeros(layer.weights.dim()),
                    scale: Array1::zeros(layer.scale.len()),
                    shift: Array1::zeros(layer.shift.len()),
                })
                .collect(),
            last: Array2::zeros(self.last.dim()),
            last_bias: Array1::zeros(WIDTH),
        }
    }

    pub fn slices(&self) -> Vec<&[f32]> {
        let mut out = Vec::new();
        for layer in &self.layers {
            out.extend([flat(&layer.weights), flat(&layer.scale), flat(&layer.shift)]);
        }
        out.extend([flat(&self.last), flat(&self.last_bias)]);
        out
    }

    pub fn slices_mut(&mut self) -> Vec<&mut [f32]> {
        let mut out = Vec::new();
        for layer in &mut self.layers {
            out.extend([
                flat_mut(&mut layer.weights),
                flat_mut(&mut layer.scale),
                flat_mut(&mut layer.shift),
            ]);
        }
        out.extend([flat_mut(&mut self.last), flat_mut(&mut self.last_bias)]);
        out
    }
}

fn flat<D: ndarray::Dimension>(values: &ndarray::Array<f32, D>) -> &[f32] {
    values.as_slice().expect("parameters are contiguous")
}

fn flat_mut<D: ndarray::Dimension>(values: &mut ndarray::Array<f32, D>) -> &mut [f32] {
    values.as_slice_mut().expect("parameters are contiguous")
}

pub struct Net {
    pub params: Params,
    running: Vec<Running>,
}

pub struct Pass {
    input: Array2<f32>,
    steps: Vec<Normed>,
    pub output: Array2<f32>,
    lengths: Vec<f32>,
    kept: f32,
}

impl Pass {
    fn below(&self, at: usize) -> ArrayView2<'_, f32> {
        match at {
            0 => self.input.view(),
            _ => self.steps[at - 1].active.view(),
        }
    }
}

impl Net {
    pub fn new(input: usize, rng: &mut StdRng) -> Self {
        let running = LAYERS.iter().map(|&width| Running::new(width)).collect();
        Self {
            params: Params::new(input, rng),
            running,
        }
    }

    pub fn forward(&mut self, input: Array2<f32>) -> Pass {
        self.forward_dropping(input, None)
    }

    pub fn forward_dropping(&mut self, input: Array2<f32>, drop: Option<(f32, u64)>) -> Pass {
        let rows = input.nrows();
        let mut steps: Vec<Normed> = Vec::with_capacity(LAYERS.len());
        for (at, (layer, running)) in self.params.layers.iter().zip(&mut self.running).enumerate() {
            let below = steps.last().map_or(input.view(), |step| step.active.view());
            let raw = fresh(1.0, below, layer.weights.t());
            let drop = drop.filter(|_| at == 0);
            steps.push(norm::normalise(
                raw,
                running,
                (&layer.scale, &layer.shift),
                drop,
            ));
        }
        let top = steps.last().map_or(input.view(), |step| step.active.view());
        let mut output = Array2::from_shape_fn((rows, WIDTH), |(_, at)| self.params.last_bias[at]);
        product(1.0, top, self.params.last.t(), 1.0, output.view_mut());
        let lengths = norm::unit(&mut output);
        Pass {
            input,
            steps,
            output,
            lengths,
            kept: norm::dropout_scale(drop),
        }
    }

    pub fn backward(&self, pass: &Pass, gradient: Array2<f32>) -> Params {
        self.backward_parts(pass, gradient, false).0
    }

    pub fn backward_input(&self, pass: &Pass, gradient: Array2<f32>) -> (Params, Array2<f32>) {
        let (grads, input) = self.backward_parts(pass, gradient, true);
        (grads, input.expect("asked for the input gradient"))
    }

    fn backward_parts(
        &self,
        pass: &Pass,
        mut gradient: Array2<f32>,
        want_input: bool,
    ) -> (Params, Option<Array2<f32>>) {
        norm::unit_backward(&mut gradient, &pass.output, &pass.lengths);
        let layers = self.params.layers.len();
        let last = later(1.0, gradient.t(), pass.below(layers));
        let last_bias = norm::sums(&gradient, None);
        let mut below = fresh(1.0, gradient.view(), self.params.last.view());
        let mut grads = Vec::with_capacity(layers);
        let mut input = None;
        for at in (0..layers).rev() {
            let (layer, step) = (&self.params.layers[at], &pass.steps[at]);
            let kept = if at == 0 { pass.kept } else { 1.0 };
            let (scale, shift) = norm::backward(&mut below, step, &layer.scale, kept);
            let next =
                (at > 0 || want_input).then(|| fresh(1.0, below.view(), layer.weights.view()));
            grads.push((later_turned(1.0, below, pass.below(at)), scale, shift));
            match next {
                Some(next) if at > 0 => below = next,
                next => {
                    input = next;
                    break;
                }
            }
        }
        let layers = grads
            .into_iter()
            .rev()
            .map(|(weights, scale, shift)| Layer {
                weights: weights.wait(),
                scale,
                shift,
            })
            .collect();
        let grads = Params {
            layers,
            last: last.wait(),
            last_bias,
        };
        (grads, input)
    }

    pub fn embed(&self, input: ArrayView2<f32>) -> Array2<f32> {
        let mut hidden: Option<Array2<f32>> = None;
        for (layer, running) in self.params.layers.iter().zip(&self.running) {
            let below = hidden.as_ref().map_or(input, |table| table.view());
            let raw = fresh(1.0, below, layer.weights.t());
            hidden = Some(norm::infer(raw, running, (&layer.scale, &layer.shift)));
        }
        let top = hidden.as_ref().map_or(input, |table| table.view());
        let mut output =
            Array2::from_shape_fn((input.nrows(), WIDTH), |(_, at)| self.params.last_bias[at]);
        product(1.0, top, self.params.last.t(), 1.0, output.view_mut());
        norm::unit(&mut output);
        output
    }
}
