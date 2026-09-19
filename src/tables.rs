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
                "the composition table holds {} contigs and the coverage table {n_contigs}",
                tnf.kmer_table.nrows()
            );
        }

        debug!("Filtering TNF table.");
        tnf.filter_by_name(&filtered)?;
        if tnf.kmer_table.nrows() != coverage.table.nrows() {
            bail!("the two tables hold different contigs after the length filter");
        }
        tnf.clr(&coverage.contig_lengths)?;
        if let Some(target) = sources.distance.kmer_pca {
            let _timer = crate::timing::scope("kmer_pca");
            crate::kmers::pca::project(&mut tnf.kmer_table, target)?;
        }

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
