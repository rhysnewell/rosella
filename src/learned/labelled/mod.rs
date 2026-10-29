pub mod corpus;
pub mod loss;
mod pieces;

use anyhow::Result;
use log::info;
use ndarray::{Array2, ArrayView2};
use rand::{SeedableRng, rngs::StdRng};

use crate::learned::corpus::{Scales, depth_of, embedded};
use crate::learned::device;
use crate::learned::draw::Draw;
use crate::learned::net::Net;
use crate::learned::optim::{AdamW, Epochs, Plateau, Run, Schedule};
use crate::progress::{self, Stage};
use corpus::{Corpus, Inputs};

const BIN_BATCH: usize = 4096;
const CROP_BATCH: usize = 2048;
const TEMPERATURE: f32 = 0.1;
const PEAK_RATE: f64 = 2e-3;
const DECAY: f64 = 1e-4;
const LAYERS: [usize; 2] = [512, 512];
const WIDTH: usize = 64;

enum Engine {
    Host(Box<AdamW>),
    Device(Box<device::Single>),
}

pub struct Labelled {
    net: Net,
    scales: Scales,
    held: Vec<usize>,
}

impl Labelled {
    pub fn train(inputs: &Inputs) -> Result<Option<Self>> {
        let corpus = Corpus::new(inputs)?;
        if corpus.labelled_bins() < 2 {
            info!(
                "Learned space: off ({} near complete bins)",
                corpus.labelled_bins()
            );
            return Ok(None);
        }
        let _timer = crate::timing::scope("learned_train");
        let mut rng = StdRng::seed_from_u64(inputs.seed);
        let mut net = Net::new(corpus.scales.width(), &LAYERS, WIDTH, &mut rng);
        let per_bin = (BIN_BATCH / corpus.labelled_bins()).max(2);
        let per_epoch = corpus.others.len().div_ceil(CROP_BATCH).max(1);
        let epochs = inputs.epochs;
        let mut schedule = Schedule::new(epochs * per_epoch, PEAK_RATE);
        let rows = corpus.bin_rows(per_bin) + 2 * CROP_BATCH;
        let mut engine = match device::Single::new(&net, rows, DECAY)? {
            Some(single) => Engine::Device(Box::new(single)),
            None => Engine::Host(Box::new(AdamW::new(&net.params, DECAY))),
        };
        let progress = progress::stepped(Stage::LearningSpace, (epochs * per_epoch) as u64);
        let lengths = corpus.others.iter().map(|row| inputs.lengths[*row]);
        let mut draw = Draw::new(&lengths.collect::<Vec<_>>(), inputs.fast);
        let offset = corpus.labelled_bins() as u32;
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
            batch: CROP_BATCH,
            offset,
            pairs: None,
            extra: &[],
        };
        let trained = shape.run(run, |rng, others, (rate, momentum)| {
            let batch = corpus.batch(rng, others, per_bin);
            let (loss, rows) = match &mut engine {
                Engine::Host(optimiser) => {
                    let pass = net.forward(batch.features);
                    let (loss, gradient, rows) =
                        loss::contrastive(pass.output.view(), &batch.labels, TEMPERATURE);
                    let grads = net.backward(&pass, gradient);
                    optimiser.update(&mut net.params, &grads, rate, momentum);
                    (loss, rows)
                }
                Engine::Device(single) => single.step(
                    &batch.features,
                    &batch.labels,
                    (rate, momentum),
                    TEMPERATURE,
                )?,
            };
            Ok(Some((loss, batch.labels, rows)))
        })?;
        if let Engine::Device(single) = engine {
            single.copy_into(&mut net)?;
        }
        progress.set_message(format!(
            "{} contigs in {} bins and {} unbinned, {} of {epochs} epochs{}, loss {:.4}",
            corpus.held_rows,
            corpus.labelled_bins(),
            corpus.others.len(),
            trained.epochs,
            draw.note(trained.drawn, trained.epochs * corpus.others.len()),
            trained.last.unwrap_or(f64::NAN)
        ));
        progress::done(&progress);
        let mut held = corpus.held;
        held.sort_unstable();
        Ok(Some(Self {
            net,
            scales: corpus.scales,
            held,
        }))
    }

    // The weight fit draws only label contigs, whose own halves were never a training pair.
    pub fn holds(&self, row: usize) -> bool {
        self.held.binary_search(&row).is_ok()
    }

    pub fn embed(&self, composition: ArrayView2<f64>, coverage: ArrayView2<f64>) -> Array2<f64> {
        self.embed_with(composition, |row| depth_of(coverage, row))
    }

    pub(crate) fn embed_with(
        &self,
        composition: ArrayView2<f64>,
        depth: impl Fn(usize) -> Vec<f64> + Sync,
    ) -> Array2<f64> {
        embedded(
            &self.scales,
            self.net.width(),
            composition,
            depth,
            |input| self.net.embed(input),
        )
    }
}
