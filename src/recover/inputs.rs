use std::path;

use anyhow::Result;
use log::debug;

use crate::{
    cli::RecoverArgs, clustering::graph_partition::Partition,
    coverage::coverage_table::CoverageTable, embedding::metrics::DistanceSettings,
    kmers::kmer_counting::KmerFrequencyTable, kmers::sketch::ContigSketches,
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
}

/// The search runs between the other stages rather than beside them. Overlapping it with
/// coverage bought wall by asking for more threads than the box has.
fn run_search(
    args: &RecoverArgs,
    assembly: &str,
    floor: usize,
) -> Result<crate::markers::MarkerAnnotation> {
    let rules = crate::markers::MarkerRules {
        fragment_span: args.markers.marker_fragment_span,
    };
    let built = crate::markers::MarkerAnnotation::build(
        assembly,
        floor,
        args.runtime.threads,
        args.markers.hmm_shards.map(usize::from),
        rules,
        args.markers.marker_cache.as_deref().map(path::Path::new),
    )?;
    if let Some(path) = &args.reports.marker_report {
        built.report(path::Path::new(path))?;
    }
    Ok(built)
}

pub fn read_inputs(args: &RecoverArgs) -> Result<Inputs> {
    let output_directory = args.common.output_directory.clone();
    let assembly = args.assembly.clone();
    let cutoff = args.binning.min_contig_size;
    let min_contig_size = args
        .binning
        .attractor_floor
        .map_or(cutoff, |floor| floor.min(cutoff));
    let tables = crate::tables::Tables::build(&crate::tables::Sources {
        assembly: &assembly,
        common: &args.common,
        min_contig_size,
        coverage: &args.coverage,
        mapping: &args.mapping,
        filtering: &args.filtering,
        alignment: &args.alignment,
        trimming: &args.trimming,
        distance: &args.distance,
        threads: args.runtime.threads,
    })?;
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
    let quality = run_search(args, &assembly, min_contig_size)?
        .select(&coverage_table.contig_names)?
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
    })
}

// Rows at or above the cutoff come first and in assembly order, so a pass over them alone sees
// the same rows, and so the same seeds, as a run with no shorter contigs.
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
