use std::path;

use anyhow::Result;
use log::{debug, info};

use crate::{
    cli::{Length, RecoverArgs},
    markers::{Annotator, ContigMarkers},
    recover::settings::Settings,
    tables::{Lengths, Sources, Tables},
};

pub struct Inputs {
    pub tables: Tables,
    pub annotator: Annotator,
    pub quality: ContigMarkers,
    pub oracle: Vec<Vec<usize>>,
}

pub fn sources(args: &RecoverArgs, lengths: Lengths) -> Sources<'_> {
    Sources {
        assembly: &args.assembly,
        common: &args.common,
        reads: &args.reads,
        distance: &args.distance,
        lengths,
        threads: args.runtime.threads,
    }
}

pub fn read_inputs(args: &RecoverArgs, settings: &Settings) -> Result<Inputs> {
    let cutoff = settings.lengths.cutoff;
    let given = settings.lengths.attach_floor;
    info!(
        "Partitioning contigs from {cutoff} bp{}.",
        chosen(args.binning.min_contig_size)
    );
    match given {
        Some(floor) => info!("Attaching contigs from {floor} bp, as given."),
        None => info!("Attaching shorter contigs down to where their markers turn foreign."),
    }
    // The weight halves contigs long enough that each half could be binned on its own.
    let mut tables = Tables::build(&sources(
        args,
        Lengths {
            sketch: settings.cycle.dissolve.then_some(cutoff),
            halves: Some(2 * cutoff),
            ..Lengths::from(given.unwrap_or(0), given.unwrap_or(cutoff))
        },
    ))?;
    long_first(&mut tables, cutoff, given.is_none());
    if let Some(sketches) = tables.sketches.as_mut() {
        sketches.align_to(&tables.coverage.contig_names)?;
    }
    tables.read_links(args.graph.assembly_graph.as_deref())?;

    let annotator = Annotator {
        assembly: args.assembly.clone(),
        threads: args.runtime.threads,
        shards: args.markers.hmm_shards.map(usize::from),
        rules: crate::markers::MarkerRules {
            fragment_span: args.markers.marker_fragment_span,
        },
        cache: args.markers.marker_cache.as_ref().map(path::PathBuf::from),
        checkm: args.markers.checkm,
    };
    let quality = annotator
        .annotate(given.unwrap_or(cutoff)..usize::MAX)?
        .select_present(&tables.coverage.contig_names)?
        .with_lengths(tables.coverage.contig_lengths.clone());
    let oracle = match &args.reports.dissolve_oracle {
        Some(path) => {
            let groups = crate::refine::oracle::read_groups(path, &tables.coverage.contig_names)?;
            debug!("Offering the pool {} groups from {path}.", groups.len());
            groups
        }
        None => Vec::new(),
    };

    Ok(Inputs {
        tables,
        annotator,
        quality,
        oracle,
    })
}

fn chosen(length: Length) -> &'static str {
    match length {
        Length::Auto => ", the default",
        Length::Given(_) => ", as given",
    }
}

// Long rows keep assembly order, so a pass over them alone matches a run without short contigs.
// Deferred rows are laid out longest first, as the walk takes them, and appended band by band.
fn long_first(tables: &mut Tables, cutoff: usize, defer: bool) {
    let Tables {
        coverage,
        tnf: composition,
        ..
    } = tables;
    let lengths = &coverage.contig_lengths;
    let mut order = (0..lengths.len())
        .filter(|row| lengths[*row] >= cutoff)
        .collect::<Vec<_>>();
    let loaded = if defer { order.len() } else { lengths.len() };
    let start = order.len();
    order.extend((0..lengths.len()).filter(|row| lengths[*row] < cutoff));
    if defer {
        order[start..].sort_by_key(|row| (std::cmp::Reverse(lengths[*row]), *row));
    }
    if loaded == lengths.len() && order.iter().enumerate().all(|(at, row)| at == *row) {
        return;
    }
    let rows = &order[..loaded];
    coverage.table = coverage.table.select(ndarray::Axis(0), rows);
    coverage.average_depths = crate::rows::reorder(&coverage.average_depths, rows);
    coverage.contig_names = crate::rows::reorder(&coverage.contig_names, &order);
    coverage.contig_lengths = crate::rows::reorder(&coverage.contig_lengths, &order);
    if !defer {
        composition.kmer_table = composition.kmer_table.select(ndarray::Axis(0), rows);
        composition.contig_names = crate::rows::reorder(&composition.contig_names, rows);
    }
}
