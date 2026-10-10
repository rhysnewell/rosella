use anyhow::Result;
use log::info;
use ndarray::s;

use crate::learned::{Learned, corpus::Inputs, windows::Windows};
use crate::recover::settings::Settings;
use crate::tables::Tables;

// The space reads no labels, so it trains while the markers are still being annotated. Its own
// pool keeps its row work from queueing behind gene calling in the shared one.
pub(super) fn learn_space(
    tables: &mut Tables,
    settings: &Settings,
    contigs: &[usize],
) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(rayon::current_num_threads())
        .build()?
        .install(|| learn(tables, settings, contigs))
}

fn learn(tables: &mut Tables, settings: &Settings, contigs: &[usize]) -> Result<()> {
    let lengths = &tables.coverage.contig_lengths;
    let long = contigs.partition_point(|contig| lengths[*contig] >= settings.lengths.cutoff);
    let spans = match &settings.window_depth {
        Some(dir) => {
            let _timer = crate::timing::scope("learned_windows");
            let windows = Windows::read(dir, &tables.coverage.sample_names)?;
            Some(windows.spans(
                &tables.coverage.contig_names,
                lengths,
                long,
                settings.seeds.seed,
            )?)
        }
        None => None,
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
        spans: spans.as_ref(),
    })?;
    let Some(learned) = learned else {
        return Ok(());
    };
    let rows = {
        let _timer = crate::timing::scope("learned_embed");
        learned.embed(tables.tnf.kmer_table.view(), tables.coverage.table.view())
    };
    info!(
        "Composition now reads {} learned dimensions for {} contigs.",
        rows.ncols(),
        rows.nrows()
    );
    tables.tnf.kmer_table = rows;
    tables.learned = Some(learned);
    Ok(())
}
