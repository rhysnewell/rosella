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
    homology::{Homology, homology_settings},
    kmers::kmer_counting::{KmerFrequencyTable, count_kmers},
    kmers::sketch::{ContigSketches, SketchParams},
    external::hmmer_engine::HmmerEngine,
    markers::ContigMarkers,
    recover::recover_engine::RECOVER_FASTA_EXTENSION,
    recover::settings::{distance_settings, transform_table},
};

pub struct Inputs {
    pub output_directory: String,
    pub assembly: String,
    pub min_contig_size: usize,
    pub coverage_table: CoverageTable,
    pub tnf_table: KmerFrequencyTable,
    pub sketches: Option<ContigSketches>,
    pub homology: Option<Homology>,
    pub markers: Option<ContigMarkers>,
    pub distance: DistanceSettings,
    pub partition: Partition,
}

pub fn read_inputs(args: &RecoverArgs) -> Result<Inputs> {
    // Read before the coverage stage, so a typo costs a message rather than a full run.
    let distance = distance_settings(&args.distance)?;
    if args.markers {
        HmmerEngine::check_installed()?;
    }
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
    transform_table(
        &mut tnf_table,
        distance.composition,
        &coverage_table.contig_lengths,
    )?;

    info!(
        "{} valid contigs, {} filtered contigs.",
        coverage_table.table.nrows(),
        filtered_contigs.len()
    );
    let partition = Partition::parse(&args.binning.partition)
        .expect("clap restricts the value")
        .resolve(&coverage_table.contig_lengths);
    let sketches = (!args.no_eject_duplicated || !args.no_rescue)
        .then(|| {
            let _timer = crate::timing::scope("sketch");
            info!("Sketching contig k-mers.");
            ContigSketches::build(
                &assembly,
                SketchParams {
                    kmer_size: args.duplication_kmer_size,
                    scale: args.duplication_scale,
                },
            )
            .map(|mut built| {
                built.filter_by_name(&filtered_contigs);
                built
            })
        })
        .transpose()?;
    if let Some(built) = &sketches {
        assert_eq!(
            coverage_table.table.nrows(),
            built.len(),
            "Coverage table and sketch table have different number of contigs."
        );
    }
    let homology = homology_settings(&args.binning, min_contig_size)
        .map(|settings| {
            Homology::build(
                &assembly,
                args.common.threads,
                settings,
                &coverage_table.contig_names,
                &coverage_table.contig_lengths,
            )
        })
        .transpose()?;
    let markers = args.markers
        .then(|| {
            let _timer = crate::timing::scope("markers");
            info!("Finding single copy markers.");
            ContigMarkers::annotate(
                &assembly,
                &coverage_table.contig_names,
                args.common.threads,
            )
        })
        .transpose()?;

    Ok(Inputs {
        output_directory,
        assembly,
        min_contig_size,
        coverage_table,
        tnf_table,
        sketches,
        homology,
        markers,
        distance,
        partition,
    })
}
