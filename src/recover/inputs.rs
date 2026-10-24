use std::collections::HashSet;
use std::path;
use std::thread::JoinHandle;

use anyhow::Result;
use log::{debug, info};

use crate::{
    cli::RunArgs,
    markers::{Annotator, ContigMarkers},
    recover::origin::{BinSet, Origin, Placed},
    recover::settings::Settings,
    tables::{Lengths, Sources, Tables},
};

pub struct Inputs {
    pub tables: Tables,
    pub annotator: Annotator,
    pub marking: JoinHandle<Result<ContigMarkers>>,
    pub oracle: Vec<Vec<usize>>,
    pub given: Option<Placed>,
    pub candidates: Vec<BinSet>,
}

pub fn sources<'a>(
    args: &'a RunArgs,
    lengths: Lengths,
    only: Option<&'a HashSet<String>>,
) -> Sources<'a> {
    Sources {
        assembly: &args.assembly,
        common: &args.common,
        reads: &args.reads,
        distance: &args.distance,
        lengths,
        threads: args.runtime.threads,
        only,
    }
}

pub fn read_inputs(args: &RunArgs, settings: &Settings, origin: &Origin) -> Result<Inputs> {
    let cutoff = settings.lengths.cutoff;
    let given = settings.lengths.attach_floor;
    info!("Minimum contig size: {cutoff}");
    match given {
        Some(floor) => debug!("Attach floor: {floor} bp"),
        None => debug!("Attach floor: auto"),
    }
    let universe = origin.universe();
    // The weight halves contigs long enough that each half could be binned on its own, and
    // publish sums the halves of every contig its k-mer test measures.
    let mut tables = laid_out(
        args,
        settings,
        Lengths {
            sketch: settings.cycle.dissolve.then_some(cutoff),
            halves: Some((2 * cutoff).min(crate::markers::replicon::LEAST_BASES)),
            ..row_lengths(settings)
        },
        universe.as_ref(),
    )?;
    if let Some(sketches) = tables.sketches.as_mut() {
        sketches.align_to(&tables.coverage.contig_names)?;
    }
    tables.read_links(args.graph.assembly_graph.as_deref())?;

    let annotator = Annotator {
        assembly: args.assembly.clone(),
        rules: crate::markers::MarkerRules {
            fragment_span: args.markers.marker_fragment_span,
        },
        cache: args.markers.marker_cache.as_ref().map(path::PathBuf::from),
        content: tables.content,
        checkm: args.markers.checkm,
        proteins: None,
    };
    let marking = {
        let annotator = annotator.clone();
        let names = tables.coverage.contig_names.clone();
        let lengths = tables.coverage.contig_lengths.clone();
        let floor = given.unwrap_or(cutoff);
        // A cache entry serves the whole assembly, so only a run without one calls genes on
        // the given bins alone.
        let subset = (universe.is_some() && annotator.cache.is_none()).then(|| {
            names
                .iter()
                .zip(&lengths)
                .filter(|(_, length)| **length >= floor)
                .map(|(name, _)| name.clone())
                .collect::<Vec<_>>()
        });
        std::thread::spawn(move || {
            crate::pool::install(|| {
                let annotation = match &subset {
                    Some(held) => annotator.annotate_named(held, "given.fna")?,
                    None => annotator.annotate(floor..usize::MAX)?,
                };
                Ok(annotation.select_present(&names)?.with_lengths(lengths))
            })
        })
    };
    let oracle = match &args.reports.dissolve_oracle {
        Some(path) => {
            let groups = crate::refine::oracle::read_groups(path, &tables.coverage.contig_names)?;
            debug!("Oracle groups: {} from {path}", groups.len());
            groups
        }
        None => Vec::new(),
    };
    let given = origin.place(&tables.coverage.contig_names, universe);
    let candidates = origin.candidates(&tables.coverage.contig_names);

    Ok(Inputs {
        tables,
        annotator,
        marking,
        oracle,
        given,
        candidates,
    })
}

pub fn row_lengths(settings: &Settings) -> Lengths {
    let floor = settings.lengths.attach_floor;
    Lengths::from(floor.unwrap_or(0), floor.unwrap_or(settings.lengths.cutoff))
}

pub fn laid_out(
    args: &RunArgs,
    settings: &Settings,
    lengths: Lengths,
    only: Option<&HashSet<String>>,
) -> Result<Tables> {
    let mut tables = Tables::build(&sources(args, lengths, only))?;
    let defer = settings.lengths.attach_floor.is_none();
    long_first(&mut tables, settings.lengths.cutoff, defer);
    Ok(tables)
}

// Rows are laid out by a hash of the contig name, so every seed and tie-break that reads a row
// number lands on the same contig whatever order the assembly listed them in. Deferred rows go
// longest first, as the walk takes them.
pub fn row_order(names: &[String], lengths: &[usize], cutoff: usize, defer: bool) -> Vec<usize> {
    let keys = names
        .iter()
        .map(|name| crate::digest::mix(crate::digest::fold(name.as_bytes())))
        .collect::<Vec<_>>();
    let key = |row: &usize| (keys[*row], &names[*row]);
    let (mut order, mut short): (Vec<usize>, Vec<usize>) =
        (0..lengths.len()).partition(|row| lengths[*row] >= cutoff);
    order.sort_by_key(key);
    match defer {
        true => short.sort_by_key(|row| (std::cmp::Reverse(lengths[*row]), key(row))),
        false => short.sort_by_key(key),
    }
    order.extend(short);
    order
}

fn long_first(tables: &mut Tables, cutoff: usize, defer: bool) {
    let Tables {
        coverage,
        tnf: composition,
        ..
    } = tables;
    let lengths = &coverage.contig_lengths;
    let order = row_order(&coverage.contig_names, lengths, cutoff, defer);
    let long = lengths.iter().filter(|length| **length >= cutoff).count();
    let loaded = if defer { long } else { lengths.len() };
    if loaded == lengths.len() && order.iter().enumerate().all(|(at, row)| at == *row) {
        return;
    }
    // Deferred, the composition holds only the long rows, in table order.
    let composed = match defer {
        true => {
            let mut slot = vec![0; lengths.len()];
            for (at, row) in (0..lengths.len())
                .filter(|row| lengths[*row] >= cutoff)
                .enumerate()
            {
                slot[row] = at;
            }
            order[..loaded]
                .iter()
                .map(|row| slot[*row])
                .collect::<Vec<_>>()
        }
        false => order.clone(),
    };
    let rows = &order[..loaded];
    coverage.table = coverage.table.select(ndarray::Axis(0), rows);
    coverage.average_depths = crate::rows::reorder(&coverage.average_depths, rows);
    coverage.contig_names = crate::rows::reorder(&coverage.contig_names, &order);
    coverage.contig_lengths = crate::rows::reorder(&coverage.contig_lengths, &order);
    composition.take_rows(&composed);
}
