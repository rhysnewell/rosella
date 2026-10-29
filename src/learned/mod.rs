pub mod corpus;
pub mod device;
pub mod draw;
pub mod gemm;
pub mod labelled;
pub mod loss;
mod neighbours;
pub mod net;
mod norm;
pub mod optim;
pub mod pieces;
mod sample;
pub mod store;
pub mod towers;
pub mod windows;

use anyhow::Result;
use log::info;
use ndarray::{Array2, ArrayView2};
use rand::{Rng, SeedableRng, rngs::StdRng};

use crate::progress::{self, Stage};
use corpus::{Corpus, Inputs, Scales, depth_of, embedded};
use draw::Draw;
use labelled::{Labelled, corpus::half_depths};
use optim::{Epochs, Plateau, Run, Schedule};
use towers::Towers;

pub(crate) const VIEWS_PER_STEP: usize = 4096;
// Each positive pair against the negatives alone at 0.07, with dropout 0.2, beat the supervised
// loss at 0.1 (Khosla et al. 2020) on every tier of mouse gut.
const TEMPERATURE: f32 = 0.07;
const DROPOUT: f32 = 0.2;
const PEAK_RATE: f64 = 2e-3;
const DECAY: f64 = 1e-4;

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

    pub fn halves(
        &self,
        [first, second]: [Array2<f64>; 2],
        coverage: &[&[f64]],
        lengths: &[usize],
        seed: u64,
    ) -> [Array2<f64>; 2] {
        let [early, late] = [0, 1].map(|side| half_depths(coverage, lengths, seed, side));
        let embed = |composition: Array2<f64>, depth: Vec<Vec<f64>>| {
            let depth = move |row: usize| depth[row].clone();
            match self {
                Self::Windowed(learned) => learned.embed_with(composition.view(), depth),
                Self::Labelled(labelled) => labelled.embed_with(composition.view(), depth),
            }
        };
        [embed(first, early), embed(second, late)]
    }
}

impl Learned {
    pub fn new(net: Towers, scales: corpus::Scales) -> Self {
        Self { net, scales }
    }

    pub fn train(inputs: &Inputs) -> Result<Option<Self>> {
        let stretches = inputs.stretches;
        if inputs.composition.nrows() < 1 + stretches {
            info!(
                "Learned space: off ({} contigs)",
                inputs.composition.nrows()
            );
            return Ok(None);
        }
        let per_step = VIEWS_PER_STEP / (1 + stretches);
        let per_epoch = inputs.composition.nrows().div_ceil(per_step);
        let epochs = inputs.epochs;
        let progress = progress::stepped(Stage::LearningSpace, (epochs * per_epoch) as u64);
        progress.set_message("reading contigs");
        let corpus = Corpus::new(inputs)?;
        let _timer = crate::timing::scope("learned_train");
        let mut rng = StdRng::seed_from_u64(inputs.seed);
        let mut net = Towers::new(corpus.scales.width(), corpus.scales.split(), &mut rng);
        let mut schedule = Schedule::new(epochs * per_epoch, PEAK_RATE);
        let mut engine = net.engine(DECAY)?;
        let mut draw = Draw::new(inputs.lengths, inputs.fast);
        let mut plateau = Plateau::new(inputs.patience);
        let run = Run {
            progress: &progress,
            schedule: &mut schedule,
            draw: &mut draw,
            plateau: &mut plateau,
            rng: &mut rng,
        };
        let shape = Epochs {
            epochs,
            per_epoch,
            batch: per_step,
            offset: 0,
        };
        let trained = shape.run(run, |rng, contigs, rate| {
            if contigs.is_empty() {
                return Ok(None);
            }
            let batch = corpus.batch(rng, contigs, stretches);
            let drop = Some((DROPOUT, rng.random::<u64>()));
            let (loss, rows) = net.train_step(
                &mut engine,
                (batch.features, &batch.labels),
                drop,
                rate,
                TEMPERATURE,
            )?;
            Ok(Some((loss, batch.labels, rows)))
        })?;
        net.finish(engine)?;
        progress.set_message(format!(
            "{} contigs, {} of {epochs} epochs{}, loss {:.4}",
            corpus.contigs,
            trained.epochs,
            draw.note(trained.drawn, trained.epochs * corpus.contigs),
            trained.last.unwrap_or(f64::NAN)
        ));
        progress::done(&progress);
        Ok(Some(Self::new(net, corpus.scales)))
    }

    pub fn embed(&self, composition: ArrayView2<f64>, coverage: ArrayView2<f64>) -> Array2<f64> {
        self.embed_with(composition, |row| depth_of(coverage, row))
    }

    fn embed_with(
        &self,
        composition: ArrayView2<f64>,
        depth: impl Fn(usize) -> Vec<f64> + Sync,
    ) -> Array2<f64> {
        embedded(
            &self.scales,
            self.net.width(),
            composition,
            depth,
            |input| match device::embed(&self.net, input) {
                Ok(Some(embedded)) => embedded,
                Ok(None) => self.net.embed(input),
                Err(error) => {
                    log::warn!("GPU: embedding failed ({error:#}), embedding on the CPU");
                    self.net.embed(input)
                }
            },
        )
    }
}
