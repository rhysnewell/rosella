use std::path;

use anyhow::Result;
use log::{debug, info};

use crate::{
    cli::{Length, RecoverArgs},
    clustering::graph_partition::Partition,
    coverage::coverage_table::CoverageTable,
    embedding::metrics::DistanceSettings,
    kmers::kmer_counting::KmerFrequencyTable,
    kmers::sketch::ContigSketches,
};

pub struct Inputs {
    pub output_directory: String,
    pub assembly: String,
    pub min_contig_size: usize,
    pub coverage_table: CoverageTable,
    pub coverage_file: String,
    pub tnf_table: KmerFrequencyTable,
    pub sketches: Option<ContigSketches>,
    pub links: Option<Vec<crate::assembly_graph::Link>>,
    pub quality: crate::markers::ContigMarkers,
    pub oracle: Vec<Vec<usize>>,
    pub distance: DistanceSettings,
    pub partition: Partition,
    pub dissolve: bool,
    pub cutoff: usize,
    pub attach_given: bool,
    pub annotator: crate::markers::Annotator,
}

/// The search runs between the other stages rather than beside them. Overlapping it with
/// coverage bought wall by asking for more threads than the box has.
pub fn sources(
    args: &RecoverArgs,
    min_contig_size: usize,
    composition_from: usize,
) -> crate::tables::Sources<'_> {
    crate::tables::Sources {
        assembly: &args.assembly,
        common: &args.common,
        min_contig_size,
        composition_from,
        coverage: &args.coverage,
        mapping: &args.mapping,
        filtering: &args.filtering,
        alignment: &args.alignment,
        trimming: &args.trimming,
        distance: &args.distance,
        threads: args.runtime.threads,
    }
}

pub fn read_inputs(args: &RecoverArgs) -> Result<Inputs> {
    let output_directory = args.common.output_directory.clone();
    let assembly = args.assembly.clone();
    let cutoff = args.binning.cutoff();
    let given = args.attach_floor();
    let min_contig_size = given.unwrap_or(0);
    info!(
        "Partitioning contigs from {cutoff} bp{}.",
        chosen(args.binning.min_contig_size)
    );
    match given {
        Some(floor) => info!("Attaching contigs from {floor} bp, as given."),
        None => info!("Attaching shorter contigs down to where their markers turn foreign."),
    }
    let tables =
        crate::tables::Tables::build(&sources(args, min_contig_size, given.unwrap_or(cutoff)))?;
    let (mut coverage_table, coverage_file, mut tnf_table, distance) = (
        tables.coverage,
        tables.coverage_file,
        tables.tnf,
        tables.distance,
    );
    long_first(&mut coverage_table, &mut tnf_table, cutoff, given.is_none());

    let partition = Partition::parse(&args.binning.partition).expect("clap restricts the value");
    let dissolve = crate::recover::settings::dissolve(&args.rescue.dissolve);
    let sketches = dissolve
        .then(|| {
            let _timer = crate::timing::scope("sketch");
            debug!("Sketching contig k-mers.");
            ContigSketches::build(&assembly, cutoff).and_then(|mut built| {
                built.align_to(&coverage_table.contig_names)?;
                Ok(built)
            })
        })
        .transpose()?;
    let links = args
        .graph
        .assembly_graph
        .as_ref()
        .map(|path| crate::assembly_graph::read_links(path, &coverage_table.contig_names))
        .transpose()?;
    let annotator = crate::markers::Annotator {
        assembly: assembly.clone(),
        threads: args.runtime.threads,
        shards: args.markers.hmm_shards.map(usize::from),
        rules: crate::markers::MarkerRules {
            fragment_span: args.markers.marker_fragment_span,
        },
        cache: args.markers.marker_cache.as_ref().map(path::PathBuf::from),
    };
    let quality = annotator
        .annotate(given.unwrap_or(cutoff)..usize::MAX)?
        .select_present(&coverage_table.contig_names)?
        .with_lengths(coverage_table.contig_lengths.clone())
        .counting(crate::recover::settings::duplicates(
            &args.markers.marker_duplicates,
        ))
        .with_partials(crate::recover::settings::partials(
            &args.markers.marker_partials,
        ));
    let oracle = match &args.reports.dissolve_oracle {
        Some(path) => {
            let groups = crate::refine::oracle::read_groups(path, &coverage_table.contig_names)?;
            debug!("Offering the pool {} groups from {path}.", groups.len());
            groups
        }
        None => Vec::new(),
    };

    Ok(Inputs {
        output_directory,
        assembly,
        min_contig_size,
        coverage_table,
        coverage_file,
        tnf_table,
        sketches,
        links,
        quality,
        oracle,
        distance,
        partition,
        dissolve,
        cutoff,
        attach_given: given.is_some(),
        annotator,
    })
}

fn chosen(length: Length) -> &'static str {
    match length {
        Length::Auto => ", the default until the assembly sets it",
        Length::Given(_) => ", as given",
    }
}

// Long rows keep assembly order so a pass over them alone sees the seeds of a run without short
// contigs. A walk takes short contigs longest first, so deferred rows are laid out that way and
// each band it reaches is appended, while names and lengths cover every contig from the start.
fn long_first(
    coverage: &mut CoverageTable,
    composition: &mut KmerFrequencyTable,
    cutoff: usize,
    defer: bool,
) {
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
