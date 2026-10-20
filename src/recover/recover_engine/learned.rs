use std::path::Path;

use anyhow::Result;
use log::info;
use ndarray::{Array2, Axis};

use crate::defaults::TRAINING_FLOOR;
use crate::learned::labelled::{self, Labelled};
use crate::learned::{Learned, Space, corpus::Inputs, gemm, windows::Windows};
use crate::recover::recover_engine::RecoverEngine;
use crate::recover::settings::Settings;
use crate::tables::Tables;

struct Training {
    composition: Array2<f64>,
    coverage: Array2<f64>,
    names: Vec<String>,
    lengths: Vec<usize>,
}

// A given attach floor loads its rows in name hash order, so those are picked by length. Rows
// still to load sit longest first, so the training set takes a run off their front.
fn training(tables: &Tables, cutoff: usize) -> Result<Training> {
    let floor = TRAINING_FLOOR.min(cutoff);
    let lengths = &tables.coverage.contig_lengths;
    let loaded = tables.tnf.kmer_table.nrows();
    let rows = (0..loaded)
        .filter(|row| lengths[*row] >= floor)
        .collect::<Vec<_>>();
    let unread = loaded..loaded + lengths[loaded..].partition_point(|length| *length >= floor);
    let mut composition = tables.tnf.kmer_table.select(Axis(0), &rows);
    let mut coverage = tables.coverage.table.select(Axis(0), &rows);
    if !unread.is_empty() {
        let (more, depth) = tables.read_rows(unread.clone())?;
        composition.append(Axis(0), more.view())?;
        coverage.append(Axis(0), depth.view())?;
    }
    let rows = rows.into_iter().chain(unread).collect::<Vec<_>>();
    Ok(Training {
        composition,
        coverage,
        names: crate::rows::reorder(&tables.coverage.contig_names, &rows),
        lengths: crate::rows::reorder(lengths, &rows),
    })
}

// The space reads no labels, so it trains while the markers are still being annotated. Its own
// pool keeps its row work from queueing behind gene calling in the shared one.
pub(super) fn learn_space(tables: &mut Tables, settings: &Settings, dir: &Path) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(rayon::current_num_threads())
        .build()?
        .install(|| learn(tables, settings, dir))
}

fn learn(tables: &mut Tables, settings: &Settings, dir: &Path) -> Result<()> {
    let _gpu = gemm::gpu(settings.gpu);
    let training = training(tables, settings.lengths.cutoff)?;
    let spans = {
        let _timer = crate::timing::scope("learned_windows");
        let windows = Windows::read(dir, &tables.coverage.sample_names)?;
        windows.spans(
            &training.names,
            &training.lengths,
            training.names.len(),
            settings.seeds.seed,
        )?
    };
    let learned = Learned::train(&Inputs {
        assembly: &tables.assembly,
        names: &training.names,
        lengths: &training.lengths,
        composition: training.composition.view(),
        coverage: training.coverage.view(),
        kmer_size: tables.tnf.kmer_size(),
        seed: settings.seeds.seed,
        epochs: settings.epochs,
        patience: settings.patience,
        stretches: settings.stretches,
        fast: settings.fast_train,
        spans: &spans,
    })?;
    let Some(learned) = learned else {
        return Ok(());
    };
    let rows = {
        let _timer = crate::timing::scope("learned_embed");
        learned.embed(tables.tnf.kmer_table.view(), tables.coverage.table.view())
    };
    tables.tnf.kmer_table = rows;
    tables.learned = Some(Space::Windowed(Box::new(learned)));
    Ok(())
}

impl RecoverEngine {
    // A pass in k-mer space lends its near complete bins as labels, and the run then reads
    // composition from the space they teach.
    pub(super) fn learn_labelled(&mut self, contigs: &[usize]) -> Result<()> {
        let lengths = self.lengths();
        let cutoff = self.settings.lengths.cutoff;
        let long = contigs.partition_point(|contig| lengths[*contig] >= cutoff);
        let labels = self.first_labels(&contigs[..long])?;
        if self.settings.gpu {
            info!("GPU: the labelled space trains on the CPU");
        }
        let tables = &self.tables;
        let training = training(tables, cutoff)?;
        let learned = Labelled::train(&labelled::corpus::Inputs {
            assembly: &tables.assembly,
            names: &training.names,
            lengths: &training.lengths,
            composition: training.composition.view(),
            coverage: training.coverage.view(),
            labels: &labels,
            halves: &tables.halves,
            halves_floor: 2 * cutoff,
            kmer_size: tables.tnf.kmer_size(),
            seed: self.settings.seeds.seed,
            epochs: self.settings.epochs,
            patience: self.settings.patience,
            fast: self.settings.fast_train,
        })?;
        let Some(learned) = learned else {
            return Ok(());
        };
        let rows = {
            let _timer = crate::timing::scope("learned_embed");
            learned.embed(tables.tnf.kmer_table.view(), tables.coverage.table.view())
        };
        self.tables.tnf.kmer_table = rows;
        self.tables.learned = Some(Space::Labelled(Box::new(learned)));
        Ok(())
    }
}
