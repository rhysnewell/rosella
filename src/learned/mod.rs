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

use std::collections::HashSet;

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

// Families share a step, so the steps a set spends on its markers grow with its carriers and not
// with its family count.
pub fn packed(families: &[Vec<u32>], per_step: usize) -> Vec<Vec<u32>> {
    let mut steps = Vec::new();
    let (mut step, mut held) = (Vec::new(), HashSet::new());
    let chunks = families.iter().flat_map(|family| family.chunks(per_step));
    for chunk in chunks.filter(|chunk| chunk.len() > 1) {
        let fresh = chunk
            .iter()
            .filter(|contig| !held.contains(*contig))
            .count();
        if held.len() + fresh > per_step {
            steps.push(std::mem::take(&mut step));
            held.clear();
        }
        step.extend(chunk.iter().filter(|contig| held.insert(**contig)));
    }
    steps.push(step);
    steps.retain(|step| step.len() > 1);
    steps
}

// Carriers that fit in one step are a few strains among many genomes, and on CAMI III their steps
// moved the rest of the space more than they parted them.
pub fn marker_steps(families: &[Vec<u32>], per_step: usize) -> Vec<Vec<u32>> {
    match packed(families, per_step) {
        steps if steps.len() > 1 => steps,
        _ => Vec::new(),
    }
}

// A family left with one carrier has nothing to part it from, so it is dropped.
pub fn kept_carriers(families: &[Vec<u32>], keep: impl Fn(u32) -> bool) -> Vec<Vec<u32>> {
    families
        .iter()
        .map(|family| {
            family
                .iter()
                .copied()
                .filter(|carrier| keep(*carrier))
                .collect::<Vec<_>>()
        })
        .filter(|family| family.len() > 1)
        .collect()
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
        let chunks = marker_steps(&inputs.families, per_step);
        match chunks.is_empty() {
            false => info!(
                "Marker families: {} over {} steps an epoch",
                inputs.families.len(),
                chunks.len()
            ),
            true if !inputs.families.is_empty() => {
                info!("Marker steps: off, the carriers that pass fit in one step")
            }
            true => {}
        }
        let shape = Epochs {
            epochs,
            per_epoch,
            batch: per_step,
            offset: 0,
            extra: &chunks,
        };
        // Marker steps draw from a stream of their own, so the normal steps draw as they would
        // without them.
        let mut apart_rng = StdRng::seed_from_u64(crate::digest::mix(inputs.seed));
        let trained = shape.run(run, |rng, contigs, rate, apart| {
            if contigs.is_empty() {
                return Ok(None);
            }
            let rng = if apart { &mut apart_rng } else { rng };
            let batch = corpus.batch(rng, contigs, stretches, apart);
            let drop = Some((DROPOUT, rng.random::<u64>()));
            let (loss, rows) = net.train_step(
                &mut engine,
                (batch.features, &batch.labels, &batch.mates),
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
