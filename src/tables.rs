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
    pub composition_from: usize,
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
    pub coverage_file: String,
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
        let distance = crate::recover::settings::distance_settings();
        let output_directory = &sources.common.output_directory;
        crate::bins::refuse_used(output_directory)?;
        std::fs::create_dir_all(output_directory)?;

        debug!("Calculating contig coverages.");
        let (mut coverage, coverage_file) = {
            let _timer = crate::timing::scope("coverage");
            calculate_coverage(&CoverageInputs {
                assembly: sources.assembly,
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

        let counted = coverage
            .contig_lengths
            .iter()
            .zip(&coverage.contig_names)
            .filter(|(length, _)| **length >= sources.composition_from)
            .collect::<Vec<_>>();
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
                        counted.len(),
                        sources.composition_from,
                        sources.distance.kmer_size,
                        sources.distance.write_kmer_table,
                    )?
                }
            }
        };

        debug!("Filtering TNF table.");
        let held = counted
            .iter()
            .map(|(_, name)| name.as_str())
            .collect::<std::collections::HashSet<_>>();
        let extra = tnf
            .contig_names
            .iter()
            .filter(|name| !held.contains(name.as_str()))
            .cloned()
            .collect();
        tnf.filter_by_name(&extra)?;
        if tnf.kmer_table.nrows() != counted.len() {
            let seen = tnf
                .contig_names
                .iter()
                .map(String::as_str)
                .collect::<std::collections::HashSet<_>>();
            let stray = counted
                .iter()
                .find(|(_, name)| !seen.contains(name.as_str()));
            bail!(
                "the composition table does not hold the contigs of the coverage table, {}, so \
                 they were not built from {}",
                match stray {
                    Some((_, name)) => format!("starting with {name}"),
                    None => "and holds a contig twice".to_string(),
                },
                sources.assembly
            );
        }
        let lengths = counted
            .iter()
            .map(|(length, _)| **length)
            .collect::<Vec<_>>();
        tnf.clr(&lengths)?;

        info!(
            "{} valid contigs, {} filtered contigs.",
            coverage.table.nrows(),
            filtered.len()
        );
        Ok(Self {
            coverage,
            coverage_file,
            tnf,
            distance,
        })
    }
}
