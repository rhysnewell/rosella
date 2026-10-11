use std::path::Path;

use anyhow::Result;
use ndarray::s;

use crate::learned::labelled::{self, Labelled};
use crate::learned::{Learned, Space, corpus::Inputs, gemm, windows::Windows};
use crate::recover::recover_engine::RecoverEngine;
use crate::recover::settings::Settings;
use crate::tables::Tables;

// The space reads no labels, so it trains while the markers are still being annotated. Its own
// pool keeps its row work from queueing behind gene calling in the shared one.
pub(super) fn learn_space(
    tables: &mut Tables,
    settings: &Settings,
    contigs: &[usize],
    dir: &Path,
) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(rayon::current_num_threads())
        .build()?
        .install(|| learn(tables, settings, contigs, dir))
}

fn learn(tables: &mut Tables, settings: &Settings, contigs: &[usize], dir: &Path) -> Result<()> {
    let _gpu = gemm::gpu(settings.gpu);
    let lengths = &tables.coverage.contig_lengths;
    let long = contigs.partition_point(|contig| lengths[*contig] >= settings.lengths.cutoff);
    let spans = {
        let _timer = crate::timing::scope("learned_windows");
        let windows = Windows::read(dir, &tables.coverage.sample_names)?;
        windows.spans(
            &tables.coverage.contig_names,
            lengths,
            long,
            settings.seeds.seed,
        )?
    };
    let learned = Learned::train(&Inputs {
        assembly: &tables.assembly,
        names: &tables.coverage.contig_names,
        lengths,
        composition: tables.tnf.kmer_table.slice(s![..long, ..]),
        coverage: tables.coverage.table.slice(s![..long, ..]),
        kmer_size: tables.tnf.kmer_size(),
        seed: settings.seeds.seed,
        epochs: settings.epochs,
        stretches: settings.stretches,
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
        let tables = &self.tables;
        let learned = Labelled::train(&labelled::corpus::Inputs {
            assembly: &tables.assembly,
            names: &tables.coverage.contig_names,
            lengths: &tables.coverage.contig_lengths,
            composition: tables.tnf.kmer_table.slice(s![..long, ..]),
            coverage: tables.coverage.table.slice(s![..long, ..]),
            labels: &labels,
            halves: &tables.halves,
            halves_floor: 2 * cutoff,
            kmer_size: tables.tnf.kmer_size(),
            seed: self.settings.seeds.seed,
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
