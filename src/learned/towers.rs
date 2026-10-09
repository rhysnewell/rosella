use ndarray::{Array2, ArrayView2, Axis, concatenate, s};
use rand::rngs::StdRng;

use crate::learned::net::{HIDDEN, Net, Pass, WIDTH};
use crate::learned::optim::AdamW;

// Depth gets a network of its own whose unit output joins the composition, so a contig's depth
// pattern weighs the same whatever its abundance.
pub struct Towers {
    joint: Net,
    depth: Net,
    split: usize,
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

impl Towers {
    pub fn new(width: usize, split: usize, rng: &mut StdRng) -> Self {
        let depth = Net::new(width - split, HIDDEN, rng);
        Self {
            joint: Net::new(split + WIDTH, HIDDEN, rng),
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
        concatenate(Axis(1), &[input.slice(s![.., ..self.split]), depth])
            .expect("both sides hold one row per contig")
    }

    pub fn forward(&mut self, input: Array2<f32>, drop: Option<(f32, u64)>) -> TowerPass {
        let below = drop.map(|(share, key)| (share, crate::digest::mix(key ^ u64::MAX)));
        let depth = self
            .depth
            .forward_dropping(input.slice(s![.., self.split..]).to_owned(), below);
        let joined = self.joined(input.view(), depth.output.view());
        TowerPass {
            joint: self.joint.forward_dropping(joined, drop),
            depth,
        }
    }

    pub fn step(
        &mut self,
        pass: &TowerPass,
        gradient: Array2<f32>,
        optimisers: &mut Optimisers,
        (rate, momentum): (f64, f64),
    ) {
        let (grads, input) = self.joint.backward_input(&pass.joint, gradient);
        let below = self
            .depth
            .backward(&pass.depth, input.slice(s![.., self.split..]).to_owned());
        optimisers
            .joint
            .update(&mut self.joint.params, &grads, rate, momentum);
        optimisers
            .depth
            .update(&mut self.depth.params, &below, rate, momentum);
    }

    pub fn embed(&self, input: ArrayView2<f32>) -> Array2<f32> {
        let depth = self.depth.embed(input.slice(s![.., self.split..]));
        self.joint.embed(self.joined(input, depth.view()).view())
    }
}
