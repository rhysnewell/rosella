use anyhow::Result;
use ndarray::{Array2, ArrayView2, Axis, concatenate, s};
use rand::rngs::StdRng;

use crate::learned::device;
use crate::learned::loss::{Mates, per_positive_mated};
use crate::learned::net::{Net, Params, Pass};
use crate::learned::optim::AdamW;

// A 1024 first layer read level with 2048 on mouse gut and higher across 20 CAMI III sets. It
// trains a quarter faster. Fewer layers or a narrower depth tower lost bins.
const LAYERS: [usize; 4] = [1024, 1024, 512, 256];
const WIDTH: usize = 128;

// Depth gets a network of its own whose unit output joins the composition, so a contig's depth
// pattern weighs the same whatever its abundance. That output is lifted to unit scale per column,
// as each standardised k-mer column is, since at a unit row's scale it trained slower.
pub struct Towers {
    pub(crate) joint: Net,
    pub(crate) depth: Net,
    pub(crate) split: usize,
}

pub(crate) fn lift() -> f32 {
    (WIDTH as f32).sqrt()
}

// The depth tower drops its own units, so the CPU and the GPU must draw it from one key.
pub(crate) fn depth_drop(drop: Option<(f32, u64)>) -> Option<(f32, u64)> {
    drop.map(|(share, key)| (share, crate::digest::mix(key ^ u64::MAX)))
}

pub struct Optimisers {
    joint: AdamW,
    depth: AdamW,
}

pub struct TowerPass {
    joint: Pass,
    depth: Pass,
}

impl TowerPass {
    pub fn output(&self) -> ArrayView2<'_, f32> {
        self.joint.output.view()
    }
}

// The device trains its own copy of the weights and hands them back when training ends.
pub(crate) enum Engine {
    Host(Box<Optimisers>),
    Device(Box<device::Trainer>),
}

impl Towers {
    pub fn new(width: usize, split: usize, rng: &mut StdRng) -> Self {
        let depth = Net::new(width - split, &LAYERS, WIDTH, rng);
        Self {
            joint: Net::new(split + WIDTH, &LAYERS, WIDTH, rng),
            depth,
            split,
        }
    }

    pub fn optimisers(&self, decay: f64) -> Optimisers {
        Optimisers {
            joint: AdamW::new(&self.joint.params, decay),
            depth: AdamW::new(&self.depth.params, decay),
        }
    }

    fn joined(&self, input: ArrayView2<f32>, depth: ArrayView2<f32>) -> Array2<f32> {
        let mut joined = concatenate(Axis(1), &[input.slice(s![.., ..self.split]), depth])
            .expect("both sides hold one row per contig");
        joined
            .slice_mut(s![.., self.split..])
            .mapv_inplace(|value| value * lift());
        joined
    }

    pub fn forward(&mut self, input: Array2<f32>, drop: Option<(f32, u64)>) -> TowerPass {
        let below = depth_drop(drop);
        let depth = self
            .depth
            .forward_dropping(input.slice(s![.., self.split..]).to_owned(), below);
        let joined = self.joined(input.view(), depth.output.view());
        TowerPass {
            joint: self.joint.forward_dropping(joined, drop),
            depth,
        }
    }

    pub fn gradients(&self, pass: &TowerPass, gradient: Array2<f32>) -> (Params, Params) {
        let (grads, input) = self.joint.backward_input(&pass.joint, gradient, self.split);
        let lifted = input.mapv(|value| value * lift());
        (grads, self.depth.backward(&pass.depth, lifted))
    }

    pub fn step(
        &mut self,
        pass: &TowerPass,
        gradient: Array2<f32>,
        optimisers: &mut Optimisers,
        (rate, momentum): (f64, f64),
    ) {
        let (joint, depth) = self.gradients(pass, gradient);
        optimisers
            .joint
            .update(&mut self.joint.params, &joint, rate, momentum);
        optimisers
            .depth
            .update(&mut self.depth.params, &depth, rate, momentum);
    }

    pub(crate) fn engine(&self, decay: f64) -> Result<Engine> {
        if let Some(trainer) = device::Trainer::new(self, decay)? {
            return Ok(Engine::Device(Box::new(trainer)));
        }
        Ok(Engine::Host(Box::new(self.optimisers(decay))))
    }

    pub(crate) fn train_step(
        &mut self,
        engine: &mut Engine,
        (features, labels, mates): (Array2<f32>, &[u32], &Mates),
        drop: Option<(f32, u64)>,
        schedule: (f64, f64),
        temperature: f32,
    ) -> Result<(f64, Vec<(f64, usize)>)> {
        match engine {
            Engine::Host(optimisers) => {
                let pass = self.forward(features, drop);
                let (loss, gradient, rows) =
                    per_positive_mated(pass.output(), labels, mates, temperature);
                self.step(&pass, gradient, optimisers, schedule);
                Ok((loss, rows))
            }
            Engine::Device(trainer) => {
                trainer.step((&features, labels, mates), drop, schedule, temperature)
            }
        }
    }

    pub(crate) fn finish(&mut self, engine: Engine) -> Result<()> {
        match engine {
            Engine::Host(_) => Ok(()),
            Engine::Device(trainer) => trainer.copy_into(self),
        }
    }

    pub fn width(&self) -> usize {
        self.joint.width()
    }

    pub(crate) fn state(&self) -> Vec<&[f32]> {
        let mut out = self.joint.state();
        out.extend(self.depth.state());
        out
    }

    pub(crate) fn state_mut(&mut self) -> Vec<&mut [f32]> {
        let mut out = self.joint.state_mut();
        out.extend(self.depth.state_mut());
        out
    }

    pub fn embed(&self, input: ArrayView2<f32>) -> Array2<f32> {
        let depth = self.depth.embed(input.slice(s![.., self.split..]));
        self.joint.embed(self.joined(input, depth.view()).view())
    }
}
