use anyhow::{Result, bail};
use log::{debug, info};

use crate::cli::binning::DistanceParams;
use crate::cli::coverage::{
    AlignmentFlags, CoverageSource, CoverageTrimming, MappingParams, ReadFiltering,
};
use crate::cli::runtime::Common;
use crate::coverage::coverage_calculator::{CoverageInputs, calculate_coverage};
use crate::coverage::coverage_table::CoverageTable;
use crate::embedding::metrics::DistanceSettings;
use crate::kmers::kmer_counting::{KmerFrequencyTable, count_kmers};

pub struct Sources<'a> {
    pub assembly: &'a str,
    pub common: &'a Common,
    pub min_contig_size: usize,
    pub coverage: &'a CoverageSource,
    pub mapping: &'a MappingParams,
    pub filtering: &'a ReadFiltering,
    pub alignment: &'a AlignmentFlags,
    pub trimming: &'a CoverageTrimming,
    pub distance: &'a DistanceParams,
    pub threads: usize,
}

pub struct Tables {
    pub coverage: CoverageTable,
    pub tnf: KmerFrequencyTable,
    pub distance: DistanceSettings,
}

impl Tables {
    /// Both subcommands need the same row-aligned pair, so the guard, the stage timers and
    /// the alignment checks live here rather than once each and differently.
    pub fn build(sources: &Sources<'_>) -> Result<Self> {
        if !std::path::Path::new(sources.assembly).is_file() {
            bail!("no assembly file at {}", sources.assembly);
        }
        let distance = crate::recover::settings::distance_settings(sources.distance);
        let output_directory = &sources.common.output_directory;
        crate::bins::refuse_used(output_directory)?;
        std::fs::create_dir_all(output_directory)?;

        debug!("Calculating contig coverages.");
        let mut coverage = {
            let _timer = crate::timing::scope("coverage");
            calculate_coverage(&CoverageInputs {
                assembly: Some(sources.assembly),
                output_directory,
                threads: sources.threads,
                coverage: sources.coverage,
                mapping: sources.mapping,
                filtering: sources.filtering,
                alignment: sources.alignment,
                trimming: sources.trimming,
            })?
        };
        let n_contigs = coverage.table.nrows();

        let filtered = {
            let _timer = crate::timing::scope("length_filter");
            coverage.filter_by_length(sources.min_contig_size)?
        };
        if coverage.table.nrows() != n_contigs - filtered.len() {
            bail!("the length filter left the coverage table a different size than it removed");
        }
        if coverage.table.nrows() == 0 {
            bail!(
                "none of the {n_contigs} contigs in {} reach --min-contig-size {}, so there is \
                 nothing to bin",
                sources.assembly,
                sources.min_contig_size
            );
        }

        let mut tnf = {
            let _timer = crate::timing::scope("kmers");
            match &sources.common.kmer_frequency_file {
                Some(path) => {
                    debug!("Reading TNF table.");
                    KmerFrequencyTable::read(path)?
                }
                None => {
                    debug!("Calculating TNF table.");
                    count_kmers(
                        sources.assembly,
                        output_directory,
                        Some(n_contigs),
                        &sources.distance.kmer_size,
                        sources.distance.write_kmer_table,
                    )?
                }
            }
        };
        if tnf.kmer_table.nrows() != n_contigs {
            bail!(
                "the composition table holds {} contigs and the coverage table {n_contigs}, so \
                 they were not built from {}",
                tnf.kmer_table.nrows(),
                sources.assembly
            );
        }

        debug!("Filtering TNF table.");
        tnf.filter_by_name(&filtered)?;
        if tnf.kmer_table.nrows() != coverage.table.nrows() {
            let held = tnf
                .contig_names
                .iter()
                .map(String::as_str)
                .collect::<std::collections::HashSet<_>>();
            let stray = coverage
                .contig_names
                .iter()
                .find(|name| !held.contains(name.as_str()));
            bail!(
                "the two tables hold different contigs after the length filter, {}",
                match stray {
                    Some(name) => format!(
                        "starting with {name}, which the composition table \
                                           has not seen"
                    ),
                    None => "and the composition table holds the extras".to_string(),
                }
            );
        }
        tnf.clr(&coverage.contig_lengths)?;

        info!(
            "{} valid contigs, {} filtered contigs.",
            coverage.table.nrows(),
            filtered.len()
        );
        Ok(Self {
            coverage,
            tnf,
            distance,
        })
    }
}
