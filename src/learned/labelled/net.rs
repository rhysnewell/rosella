use crate::learned::gemm::product;
use crate::learned::norm::{self, Normed, Running};
use ndarray::{Array1, Array2, ArrayView2, s};
use rand::{Rng, rngs::StdRng};

pub const HIDDEN: usize = 512;
pub const WIDTH: usize = 64;

#[derive(Clone)]
pub struct Params {
    first: Array2<f32>,
    first_scale: Array1<f32>,
    first_shift: Array1<f32>,
    second: Array2<f32>,
    second_scale: Array1<f32>,
    second_shift: Array1<f32>,
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
        let first = uniform(HIDDEN, input);
        let second = uniform(HIDDEN, HIDDEN);
        let last = uniform(WIDTH, HIDDEN);
        let last_bias = uniform(1, HIDDEN).slice(s![0, ..WIDTH]).to_owned();
        Self {
            first,
            first_scale: Array1::ones(HIDDEN),
            first_shift: Array1::zeros(HIDDEN),
            second,
            second_scale: Array1::ones(HIDDEN),
            second_shift: Array1::zeros(HIDDEN),
            last,
            last_bias,
        }
    }

    pub fn zeros_like(&self) -> Self {
        Self {
            first: Array2::zeros(self.first.dim()),
            first_scale: Array1::zeros(HIDDEN),
            first_shift: Array1::zeros(HIDDEN),
            second: Array2::zeros(self.second.dim()),
            second_scale: Array1::zeros(HIDDEN),
            second_shift: Array1::zeros(HIDDEN),
            last: Array2::zeros(self.last.dim()),
            last_bias: Array1::zeros(WIDTH),
        }
    }

    pub fn slices(&self) -> [&[f32]; 8] {
        [
            flat(&self.first),
            flat(&self.first_scale),
            flat(&self.first_shift),
            flat(&self.second),
            flat(&self.second_scale),
            flat(&self.second_shift),
            flat(&self.last),
            flat(&self.last_bias),
        ]
    }

    pub fn slices_mut(&mut self) -> [&mut [f32]; 8] {
        [
            flat_mut(&mut self.first),
            flat_mut(&mut self.first_scale),
            flat_mut(&mut self.first_shift),
            flat_mut(&mut self.second),
            flat_mut(&mut self.second_scale),
            flat_mut(&mut self.second_shift),
            flat_mut(&mut self.last),
            flat_mut(&mut self.last_bias),
        ]
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
    running: [Running; 2],
}

pub struct Pass {
    input: Array2<f32>,
    steps: [Normed; 2],
    pub output: Array2<f32>,
    lengths: Vec<f32>,
}

impl Net {
    pub fn new(input: usize, rng: &mut StdRng) -> Self {
        Self {
            params: Params::new(input, rng),
            running: [Running::new(HIDDEN), Running::new(HIDDEN)],
        }
    }

    pub fn forward(&mut self, input: Array2<f32>) -> Pass {
        let rows = input.nrows();
        let mut first = Array2::zeros((rows, HIDDEN));
        product(
            1.0,
            input.view(),
            self.params.first.t(),
            0.0,
            first.view_mut(),
        );
        let first = norm::normalise(
            first,
            &mut self.running[0],
            (&self.params.first_scale, &self.params.first_shift),
            None,
        );
        let mut second = Array2::zeros((rows, HIDDEN));
        product(
            1.0,
            first.active.view(),
            self.params.second.t(),
            0.0,
            second.view_mut(),
        );
        let second = norm::normalise(
            second,
            &mut self.running[1],
            (&self.params.second_scale, &self.params.second_shift),
            None,
        );
        let mut output = Array2::from_shape_fn((rows, WIDTH), |(_, at)| self.params.last_bias[at]);
        product(
            1.0,
            second.active.view(),
            self.params.last.t(),
            1.0,
            output.view_mut(),
        );
        let lengths = norm::unit(&mut output);
        Pass {
            input,
            steps: [first, second],
            output,
            lengths,
        }
    }

    pub fn backward(&self, pass: &Pass, mut gradient: Array2<f32>) -> Params {
        let rows = gradient.nrows();
        let mut grads = self.params.zeros_like();
        norm::unit_backward(&mut gradient, &pass.output, &pass.lengths);
        product(
            1.0,
            gradient.t(),
            pass.steps[1].active.view(),
            0.0,
            grads.last.view_mut(),
        );
        grads.last_bias = norm::sums(&gradient, None);
        let mut hidden = Array2::zeros((rows, HIDDEN));
        product(
            1.0,
            gradient.view(),
            self.params.last.view(),
            0.0,
            hidden.view_mut(),
        );
        (grads.second_scale, grads.second_shift) =
            norm::backward(&mut hidden, &pass.steps[1], &self.params.second_scale, 1.0);
        product(
            1.0,
            hidden.t(),
            pass.steps[0].active.view(),
            0.0,
            grads.second.view_mut(),
        );
        let mut below = Array2::zeros((rows, HIDDEN));
        product(
            1.0,
            hidden.view(),
            self.params.second.view(),
            0.0,
            below.view_mut(),
        );
        (grads.first_scale, grads.first_shift) =
            norm::backward(&mut below, &pass.steps[0], &self.params.first_scale, 1.0);
        product(
            1.0,
            below.t(),
            pass.input.view(),
            0.0,
            grads.first.view_mut(),
        );
        grads
    }

    pub fn embed(&self, input: ArrayView2<f32>) -> Array2<f32> {
        let first = self.layer(input, &self.params.first, 0);
        let second = self.layer(first.view(), &self.params.second, 1);
        let mut output =
            Array2::from_shape_fn((input.nrows(), WIDTH), |(_, at)| self.params.last_bias[at]);
        product(
            1.0,
            second.view(),
            self.params.last.t(),
            1.0,
            output.view_mut(),
        );
        norm::unit(&mut output);
        output
    }

    fn layer(&self, input: ArrayView2<f32>, weights: &Array2<f32>, at: usize) -> Array2<f32> {
        let mut hidden = Array2::zeros((input.nrows(), HIDDEN));
        product(1.0, input, weights.t(), 0.0, hidden.view_mut());
        let (scale, shift) = match at {
            0 => (&self.params.first_scale, &self.params.first_shift),
            _ => (&self.params.second_scale, &self.params.second_shift),
        };
        norm::infer(hidden, &self.running[at], (scale, shift))
    }
}
