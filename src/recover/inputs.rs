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
fn run_search(
    args: &RecoverArgs,
    annotator: &crate::markers::Annotator,
    floor: usize,
) -> Result<crate::markers::MarkerAnnotation> {
    let built = annotator.annotate(floor..usize::MAX)?;
    if let Some(path) = &args.reports.marker_report {
        built.report(path::Path::new(path))?;
    }
    Ok(built)
}

pub fn sources(args: &RecoverArgs, min_contig_size: usize) -> crate::tables::Sources<'_> {
    crate::tables::Sources {
        assembly: &args.assembly,
        common: &args.common,
        min_contig_size,
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
    let tables = crate::tables::Tables::build(&sources(args, min_contig_size))?;
    let (mut coverage_table, mut tnf_table, distance) =
        (tables.coverage, tables.tnf, tables.distance);
    long_first(&mut coverage_table, &mut tnf_table, cutoff);

    let partition = Partition::parse(&args.binning.partition).expect("clap restricts the value");
    let dissolve = crate::recover::settings::dissolve(&args.rescue.dissolve);
    let sketches = dissolve
        .then(|| {
            let _timer = crate::timing::scope("sketch");
            debug!("Sketching contig k-mers.");
            ContigSketches::build(&assembly).and_then(|mut built| {
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
    let quality = run_search(args, &annotator, given.unwrap_or(cutoff))?
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
// contigs.
fn long_first(coverage: &mut CoverageTable, composition: &mut KmerFrequencyTable, cutoff: usize) {
    let lengths = &coverage.contig_lengths;
    let order = (0..lengths.len())
        .filter(|row| lengths[*row] >= cutoff)
        .chain((0..lengths.len()).filter(|row| lengths[*row] < cutoff))
        .collect::<Vec<_>>();
    if order.iter().enumerate().all(|(at, row)| at == *row) {
        return;
    }
    coverage.table = coverage.table.select(ndarray::Axis(0), &order);
    coverage.average_depths = crate::rows::reorder(&coverage.average_depths, &order);
    coverage.contig_names = crate::rows::reorder(&coverage.contig_names, &order);
    coverage.contig_lengths = crate::rows::reorder(&coverage.contig_lengths, &order);
    composition.kmer_table = composition.kmer_table.select(ndarray::Axis(0), &order);
    composition.contig_names = crate::rows::reorder(&composition.contig_names, &order);
}
