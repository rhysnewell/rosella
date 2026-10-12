pub mod corpus;
pub mod gemm;
pub mod labelled;
pub mod loss;
#[cfg(target_os = "macos")]
mod metal;
mod neighbours;
pub mod net;
pub mod optim;
pub mod pieces;
mod sample;
mod towers;
pub mod windows;

use anyhow::Result;
use log::{debug, info};
use ndarray::{Array2, ArrayView2, Axis};
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::progress::{self, Stage};
use corpus::{Corpus, Inputs, Scales, centred_unit, depth_of};
use labelled::Labelled;
use optim::{Plateau, Schedule};
use towers::Towers;

const VIEWS_PER_STEP: usize = 4096;
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

pub enum Space {
    Windowed(Box<Learned>),
    Labelled(Box<Labelled>),
}

impl Space {
    pub fn embed(&self, composition: ArrayView2<f64>, coverage: ArrayView2<f64>) -> Array2<f64> {
        match self {
            Self::Windowed(learned) => learned.embed(composition, coverage),
            Self::Labelled(labelled) => labelled.embed(composition, coverage),
        }
    }
}

impl Learned {
    pub fn train(inputs: &Inputs) -> Result<Option<Self>> {
        let stretches = inputs.stretches;
        if inputs.composition.nrows() < 1 + stretches {
            info!(
                "{} contigs to partition, too few to learn a space from, so composition stays \
                 k-mers.",
                inputs.composition.nrows()
            );
            return Ok(None);
        }
        let per_step = VIEWS_PER_STEP / (1 + stretches);
        let per_epoch = inputs.composition.nrows().div_ceil(per_step);
        let epochs = inputs.epochs;
        let progress = progress::stepped(Stage::LearningSpace, (epochs * per_epoch) as u64);
        progress.set_message("cutting stretches");
        let corpus = Corpus::new(inputs)?;
        let _timer = crate::timing::scope("learned_train");
        let mut rng = StdRng::seed_from_u64(inputs.seed);
        let mut net = Towers::new(corpus.scales.width(), corpus.scales.split(), &mut rng);
        let mut schedule = Schedule::new(epochs * per_epoch, PEAK_RATE);
        let mut optimisers = net.optimisers(DECAY);
        let mut order = (0..corpus.contigs as u32).collect::<Vec<_>>();
        let mut plateau = Plateau::new(inputs.patience);
        let (mut last, mut trained, mut ending) = (None, 0, false);
        for epoch in 0..epochs {
            progress.set_message(match last {
                Some(loss) => format!("e:{}/{epochs} l:{loss:.4}{}", epoch + 1, plateau.shown()),
                None => format!("e:{}/{epochs}", epoch + 1),
            });
            order.shuffle(&mut rng);
            for (step, contigs) in order.chunks(per_step).enumerate() {
                let batch = corpus.batch(&mut rng, contigs, stretches);
                let drop = Some((DROPOUT, rng.random::<u64>()));
                let pass = net.forward(batch.features, drop);
                let (loss, gradient) =
                    loss::per_positive(pass.output(), &batch.labels, TEMPERATURE);
                let at = schedule.at(epoch * per_epoch + step);
                net.step(&pass, gradient, &mut optimisers, at);
                plateau.record(loss);
                progress.inc(1);
            }
            let (loss, flat) = plateau.close(schedule.annealing(epoch * per_epoch));
            debug!("Epoch {epoch} loss {loss:.4}.");
            (last, trained) = (Some(loss), epoch + 1);
            if ending {
                break;
            }
            if flat && trained + 1 < epochs {
                schedule.anneal_from(trained * per_epoch, per_epoch);
                progress.set_length(((trained + 1) * per_epoch) as u64);
                ending = true;
            }
        }
        progress.set_message(format!(
            "{} contigs with window depth, {trained} of {epochs} epochs, loss {:.4}",
            corpus.contigs,
            last.unwrap_or(f64::NAN)
        ));
        progress::done(&progress);
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
