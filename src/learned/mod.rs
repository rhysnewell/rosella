pub mod corpus;
pub mod gemm;
pub mod loss;
mod neighbours;
pub mod net;
mod optim;
pub mod pieces;
mod sample;
mod towers;
pub mod windows;

use anyhow::Result;
use log::info;
use ndarray::{Array2, ArrayView2, Axis};
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use corpus::{Corpus, Inputs, Scales, centred_unit, depth_of};
use optim::Schedule;
use towers::Towers;

const VIEWS_PER_STEP: usize = 4096;
// Three stretches a step read inside the noise of five on mouse gut, at two thirds of the time.
const STRETCHES: usize = 3;
// Each positive pair against the negatives alone at 0.07, with dropout 0.2, beat the supervised
// loss at 0.1 (Khosla et al. 2020) on every tier of mouse gut.
const TEMPERATURE: f32 = 0.07;
const DROPOUT: f32 = 0.2;
const PEAK_RATE: f64 = 2e-3;
const DECAY: f64 = 1e-4;
const EMBED_ROWS: usize = 1 << 16;

pub struct Learned {
    net: Towers,
    scales: Scales,
}

impl Learned {
    pub fn train(inputs: &Inputs) -> Result<Option<Self>> {
        if inputs.composition.nrows() < 1 + STRETCHES {
            info!(
                "{} contigs to partition, too few to learn a space from, so composition stays \
                 k-mers.",
                inputs.composition.nrows()
            );
            return Ok(None);
        }
        let corpus = Corpus::new(inputs)?;
        let _timer = crate::timing::scope("learned_train");
        let mut rng = StdRng::seed_from_u64(inputs.seed);
        let mut net = Towers::new(corpus.scales.width(), corpus.scales.split(), &mut rng);
        let per_step = VIEWS_PER_STEP / (1 + STRETCHES);
        let per_epoch = corpus.contigs.div_ceil(per_step);
        let epochs = inputs.epochs;
        let schedule = Schedule::new(epochs * per_epoch, PEAK_RATE);
        let mut optimisers = net.optimisers(DECAY);
        info!(
            "Learning a space from {} contigs{}, {epochs} epochs of {per_epoch} steps.",
            corpus.contigs,
            match inputs.spans {
                Some(_) => " with window depth",
                None => "",
            }
        );
        let mut order = (0..corpus.contigs as u32).collect::<Vec<_>>();
        for epoch in 0..epochs {
            order.shuffle(&mut rng);
            let mut total = 0.0;
            for (step, contigs) in order.chunks(per_step).enumerate() {
                let batch = corpus.batch(&mut rng, contigs, STRETCHES);
                let drop = Some((DROPOUT, rng.random::<u64>()));
                let pass = net.forward(batch.features, drop);
                let (loss, gradient) =
                    loss::per_positive(pass.output(), &batch.labels, TEMPERATURE);
                let at = schedule.at(epoch * per_epoch + step);
                net.step(&pass, gradient, &mut optimisers, at);
                total += loss;
            }
            info!("Epoch {epoch} loss {:.4}.", total / per_epoch as f64);
        }
        Ok(Some(Self {
            net,
            scales: corpus.scales,
        }))
    }

    // Rows come back centred and unit length, so the composition distance 1 - rho reads 1 - cos.
    pub fn embed(&self, composition: ArrayView2<f64>, coverage: ArrayView2<f64>) -> Array2<f64> {
        let columns = self.scales.width();
        let mut out = Array2::zeros((composition.nrows(), net::WIDTH));
        for (at, (chunk, mut target)) in composition
            .axis_chunks_iter(Axis(0), EMBED_ROWS)
            .zip(out.axis_chunks_iter_mut(Axis(0), EMBED_ROWS))
            .enumerate()
        {
            let mut input = Array2::zeros((chunk.nrows(), columns));
            input
                .as_slice_mut()
                .expect("a fresh table is contiguous")
                .par_chunks_mut(columns)
                .enumerate()
                .for_each(|(row, values)| {
                    let depth = depth_of(coverage, at * EMBED_ROWS + row);
                    self.scales
                        .features(chunk.row(row).iter().copied(), &depth, values)
                });
            let embedded = centred_unit(self.net.embed(input.view()));
            target.assign(&embedded.mapv(f64::from));
        }
        out
    }
}
