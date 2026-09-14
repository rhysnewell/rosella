use std::path;

use anyhow::Result;
use log::{debug, info};

use crate::{
    cli::RecoverArgs,
    clustering::graph_partition::Partition,
    coverage::{
        coverage_calculator::{CoverageInputs, calculate_coverage},
        coverage_table::CoverageTable,
    },
    embedding::metrics::DistanceSettings,
    kmers::kmer_counting::{KmerFrequencyTable, count_kmers},
    recover::recover_engine::RECOVER_FASTA_EXTENSION,
    recover::settings::distance_settings,
};

pub struct Inputs {
    pub output_directory: String,
    pub assembly: String,
    pub min_contig_size: usize,
    pub coverage_table: CoverageTable,
    pub tnf_table: KmerFrequencyTable,
    pub links: Option<Vec<(usize, usize)>>,
    pub quality: crate::markers::ContigMarkers,
    pub oracle: Vec<Vec<usize>>,
    pub distance: DistanceSettings,
    pub partition: Partition,
    pub dissolve: bool,
}

/// The search runs between the other stages rather than beside them. Overlapping it with
/// coverage bought wall by asking for more threads than the box has.
fn run_search(args: &RecoverArgs, assembly: &str) -> Result<crate::markers::MarkerAnnotation> {
    let genes = crate::quality::orfs::GeneRules {
        min_length: args.gene_min_length,
        model_depth: args.gene_model_depth,
    };
    let rules = crate::markers::MarkerRules {
        fragment_span: args.marker_fragment_span,
        bar_offset: args.marker_bar_offset,
    };
    let built = crate::markers::MarkerAnnotation::build(
        assembly,
        args.binning.min_contig_size,
        genes,
        args.common.threads,
        args.hmm_shards.map(usize::from),
        rules,
        args.marker_cache.as_deref().map(path::Path::new),
    )?;
    if let Some(path) = &args.marker_report {
        built.report(path::Path::new(path))?;
    }
    Ok(built)
}

pub fn read_inputs(args: &RecoverArgs) -> Result<Inputs> {
    // Read before the coverage stage, so a typo costs a message rather than a full run.
    let distance = distance_settings(&args.distance)?;
    let output_directory = args.common.output_directory.clone();
    let output_directory_path = path::Path::new(&output_directory);
    if output_directory_path.exists() {
        let output_directory_files = output_directory_path.read_dir()?;
        for file in output_directory_files {
            let file = file?;
            let file_name = file.file_name();
            let file_name = file_name.to_str().unwrap();
            if file_name.ends_with(RECOVER_FASTA_EXTENSION) {
                return Err(anyhow::anyhow!(
                    "Output directory contains .fna files. Please remove them before running rosella recover."
                ));
            }
        }
    }

    let assembly = args.assembly.clone();
    std::fs::create_dir_all(&output_directory)?;
    info!("Calculating contig coverages.");
    let min_contig_size = args.binning.min_contig_size;
    let mut coverage_table = {
        let _timer = crate::timing::scope("coverage");
        calculate_coverage(&CoverageInputs {
            assembly: Some(&assembly),
            output_directory: &output_directory,
            threads: args.common.threads,
            coverage: &args.coverage,
            mapping: &args.mapping,
            filtering: &args.filtering,
            alignment: &args.alignment,
            trimming: &args.trimming,
        })?
    };
    let n_contigs = coverage_table.table.nrows();

    let filtered_contigs = {
        let _timer = crate::timing::scope("length_filter");
        coverage_table.filter_by_length(min_contig_size)?
    };
    if args.distance.ignore_coverage_variance {
        coverage_table.clear_variances();
    }

    assert_eq!(
        coverage_table.table.nrows(),
        n_contigs - filtered_contigs.len(),
        "Coverage table row count and total contigs minus filtered contigs do not match."
    );
    let mut tnf_table = {
        let _timer = crate::timing::scope("kmers");
        if let Some(kmer_table_path) = &args.common.kmer_frequency_file {
            info!("Reading TNF table.");
            KmerFrequencyTable::read(kmer_table_path)?
        } else {
            info!("Calculating TNF table.");
            count_kmers(
                &assembly,
                &output_directory,
                Some(n_contigs),
                args.distance.kmer_size as usize,
            )?
        }
    };
    assert_eq!(
        n_contigs,
        tnf_table.kmer_table.nrows(),
        "Coverage table row count and TNF table row count do not match."
    );
    debug!("Filtering TNF table.");
    tnf_table.filter_by_name(&filtered_contigs)?;
    assert_eq!(
        coverage_table.table.nrows(),
        tnf_table.kmer_table.nrows(),
        "Coverage table and TNF table have different number of contigs."
    );
    tnf_table.clr(&coverage_table.contig_lengths)?;

    info!(
        "{} valid contigs, {} filtered contigs.",
        coverage_table.table.nrows(),
        filtered_contigs.len()
    );
    let partition = Partition::parse(&args.binning.partition)
        .expect("clap restricts the value")
        .resolve();
    let dissolve = crate::recover::settings::dissolve(&args.dissolve);
    let links = args
        .assembly_graph
        .as_ref()
        .map(|path| crate::assembly_graph::read_links(path, &coverage_table.contig_names))
        .transpose()?;
    let quality = run_search(args, &assembly)?.select(&coverage_table.contig_names)?;
    let oracle = match &args.dissolve_oracle {
        Some(path) => {
            let groups = crate::refine::oracle::read_groups(path, &coverage_table.contig_names)?;
            info!("Offering the pool {} groups from {path}.", groups.len());
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
        links,
        quality,
        oracle,
        distance,
        partition,
        dissolve,
    })
}
