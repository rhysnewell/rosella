use anyhow::Result;
use log::{info, warn};

use crate::embedding::features::row_slice;
use crate::embedding::metrics::likelihood::{CoverageDistance, Evidence};
use crate::embedding::weight::composition::Halves;
use crate::embedding::weight::depth::{self, POOL, QUERIES};
use crate::recover::settings::Settings;
use crate::tables::Tables;

pub(super) fn evidence(
    tables: &Tables,
    settings: &Settings,
    contigs: &[usize],
) -> Result<Option<Evidence>> {
    if settings.coverage_distance == CoverageDistance::Overlap {
        return Ok(None);
    }
    let _timer = crate::timing::scope("weight");
    let lengths = &tables.coverage.contig_lengths;
    let long = contigs
        .iter()
        .copied()
        .filter(|contig| lengths[*contig] >= settings.lengths.cutoff)
        .collect::<Vec<_>>();
    let pool =
        crate::seeds::sample_positions(long.len(), POOL.min(long.len()), settings.seeds.seed)
            .into_iter()
            .map(|at| long[at])
            .collect::<Vec<_>>();
    let rows = |table| {
        pool.iter()
            .map(|contig| row_slice(table, *contig))
            .collect::<Vec<_>>()
    };
    let split = pool
        .iter()
        .copied()
        .filter(|contig| lengths[*contig] >= 2 * settings.lengths.cutoff)
        .take(QUERIES)
        .collect::<Vec<_>>();
    let names = split
        .iter()
        .map(|contig| tables.coverage.contig_names[*contig].as_str())
        .collect::<Vec<_>>();
    let [first, second] = tables.halves.composition(&names, tables.tnf.kmer_size())?;
    let halves = Halves {
        first: (0..split.len()).map(|at| row_slice(&first, at)).collect(),
        second: (0..split.len()).map(|at| row_slice(&second, at)).collect(),
        lengths: split.iter().map(|contig| lengths[*contig] as f64).collect(),
    };
    let model = depth::fit(
        &rows(&tables.coverage.table),
        &pool
            .iter()
            .map(|contig| lengths[*contig] as f64)
            .collect::<Vec<_>>(),
        &rows(&tables.tnf.kmer_table),
        &halves,
        settings.seeds.seed,
    );
    match &model {
        None => {
            warn!("Too few long contigs to fit the depth model, so depth is scored by overlap.")
        }
        Some(model) => {
            let fits = model
                .samples
                .iter()
                .map(|sample| match sample {
                    Some(sample) => format!(
                        "{:.0}:{:.3}:{}:{:.2}",
                        sample.read_length,
                        sample.spread.sqrt(),
                        sample.tail.degrees,
                        sample.absent.exp()
                    ),
                    None => "-".to_string(),
                })
                .collect::<Vec<_>>();
            info!(
                "Depth model from {} contigs, prior odds of a different genome {:.2}, read \
                     length, log depth spread, tail and absent share per sample [{}].",
                pool.len(),
                model.prior.exp(),
                fits.join(", ")
            );
        }
    }
    Ok(model)
}
