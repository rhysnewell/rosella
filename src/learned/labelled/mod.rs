pub mod corpus;
pub mod loss;
mod pieces;

use anyhow::Result;
use log::{debug, info};
use ndarray::{Array2, ArrayView2, Axis};
use rand::seq::SliceRandom;
use rand::{SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::learned::corpus::{Scales, centred_unit, depth_of};
use crate::learned::net::Net;
use crate::learned::optim::{AdamW, Plateau, Schedule, epoch_message};
use crate::progress::{self, Stage};
use corpus::{Corpus, Inputs, resampled};

const BIN_BATCH: usize = 4096;
const CROP_BATCH: usize = 2048;
const TEMPERATURE: f32 = 0.1;
const PEAK_RATE: f64 = 2e-3;
const DECAY: f64 = 1e-4;
const EMBED_ROWS: usize = 1 << 16;
const LAYERS: [usize; 2] = [512, 512];
const WIDTH: usize = 64;

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
                "Only {} bins to learn from, too few to train the learned space, so binning uses \
                 k-mer counts.",
                corpus.labelled_bins()
            );
            return Ok(None);
        }
        let _timer = crate::timing::scope("learned_train");
        let mut rng = StdRng::seed_from_u64(inputs.seed);
        let mut net = Net::new(corpus.scales.width(), &LAYERS, WIDTH, &mut rng);
        let per_bin = (BIN_BATCH / corpus.labelled_bins()).max(2);
        let per_epoch = corpus.others.div_ceil(CROP_BATCH).max(1);
        let epochs = inputs.epochs;
        let mut schedule = Schedule::new(epochs * per_epoch, PEAK_RATE);
        let mut optimiser = AdamW::new(&net.params, DECAY);
        let progress = progress::stepped(Stage::LearningSpace, (epochs * per_epoch) as u64);
        let mut order = (0..corpus.others as u32).collect::<Vec<_>>();
        let mut plateau = Plateau::new(inputs.patience);
        let (mut last, mut trained, mut ending) = (None, 0, false);
        for epoch in 0..epochs {
            progress.set_message(epoch_message(epoch, epochs, last));
            order.shuffle(&mut rng);
            for step in 0..per_epoch {
                let others = &order[(step * CROP_BATCH).min(order.len())
                    ..((step + 1) * CROP_BATCH).min(order.len())];
                let batch = corpus.batch(&mut rng, others, per_bin);
                let pass = net.forward(batch.features);
                let (loss, gradient) =
                    loss::contrastive(pass.output.view(), &batch.labels, TEMPERATURE);
                let grads = net.backward(&pass, gradient);
                let (rate, momentum) = schedule.at(epoch * per_epoch + step);
                optimiser.update(&mut net.params, &grads, rate, momentum);
                plateau.record(loss);
                progress.inc(1);
            }
            let (loss, flat) = plateau.close(schedule.annealing(epoch * per_epoch));
            debug!("Epoch {epoch} loss {loss:.4}{}.", plateau.shown());
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
            "{} contigs in {} bins and {} unbinned, {trained} of {epochs} epochs, loss {:.4}",
            corpus.bins,
            corpus.labelled_bins(),
            corpus.others,
            last.unwrap_or(f64::NAN)
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

    // Rows come back centred and unit length, so the composition distance 1 - rho reads 1 - cos.
    pub fn embed(&self, composition: ArrayView2<f64>, coverage: ArrayView2<f64>) -> Array2<f64> {
        self.embed_with(composition, |row| depth_of(coverage, row))
    }

    // Each half reads its own depth at half the span, so a pair is as noisy as two contigs that
    // long and the weight fit does not hand the space a free match.
    pub fn halves(
        &self,
        [first, second]: [Array2<f64>; 2],
        coverage: &[&[f64]],
        lengths: &[usize],
        seed: u64,
    ) -> [Array2<f64>; 2] {
        let depth = |side: u64| {
            coverage
                .iter()
                .zip(lengths)
                .enumerate()
                .map(|(at, (coverage, length))| {
                    let key = crate::digest::mix(seed ^ (2 * at as u64 + side));
                    let own = coverage.iter().step_by(2).copied().collect::<Vec<_>>();
                    resampled(&mut StdRng::seed_from_u64(key), &own, *length as f64 / 2.0)
                })
                .collect::<Vec<_>>()
        };
        let (early, late) = (depth(0), depth(1));
        [
            self.embed_with(first.view(), |row| early[row].clone()),
            self.embed_with(second.view(), |row| late[row].clone()),
        ]
    }

    fn embed_with(
        &self,
        composition: ArrayView2<f64>,
        depth: impl Fn(usize) -> Vec<f64> + Sync,
    ) -> Array2<f64> {
        let columns = self.scales.width();
        let mut out = Array2::zeros((composition.nrows(), self.net.width()));
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
                    self.scales.features(
                        chunk.row(row).iter().copied(),
                        &depth(at * EMBED_ROWS + row),
                        values,
                    )
                });
            let embedded = centred_unit(self.net.embed(input.view()));
            target.assign(&embedded.mapv(f64::from));
        }
        out
    }
}
