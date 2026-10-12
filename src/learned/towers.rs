use ndarray::{Array2, ArrayView2, Axis, concatenate, s};
use rand::rngs::StdRng;

use crate::learned::net::{Net, Pass, WIDTH};
use crate::learned::optim::AdamW;

// Depth gets a network of its own whose unit output joins the composition, so a contig's depth
// pattern weighs the same whatever its abundance. That output is lifted to unit scale per column,
// as each standardised k-mer column is, since at a unit row's scale it trained slower.
pub struct Towers {
    joint: Net,
    depth: Net,
    split: usize,
}

fn lift() -> f32 {
    (WIDTH as f32).sqrt()
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
        let depth = Net::new(width - split, rng);
        Self {
            joint: Net::new(split + WIDTH, rng),
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
        let lifted = input
            .slice(s![.., self.split..])
            .mapv(|value| value * lift());
        let below = self.depth.backward(&pass.depth, lifted);
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
