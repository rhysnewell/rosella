use std::{collections::HashSet, process::Command};

use anyhow::Result;
use log::{debug, info};

use crate::cli::{AlignmentFlags, CoverageTrimming, MappingParams, ReadFiltering};
use crate::coverage::{
    coverage_calculator::{CoverageInputs, ReadCollection},
    coverage_table::CoverageTable,
};

pub struct CovermEngine<'a> {
    assembly: &'a str,
    threads: usize,
    mapping: &'a MappingParams,
    filtering: &'a ReadFiltering,
    alignment: &'a AlignmentFlags,
    trimming: &'a CoverageTrimming,
}

impl<'a> CovermEngine<'a> {
    pub fn new(inputs: &'a CoverageInputs<'a>) -> Result<Self> {
        std::fs::create_dir_all(inputs.output_directory)?;
        check_coverm_is_installed()?;

        Ok(Self {
            assembly: inputs.assembly()?,
            threads: inputs.threads,
            mapping: inputs.mapping,
            filtering: inputs.filtering,
            alignment: inputs.alignment,
            trimming: inputs.trimming,
        })
    }

    pub fn run(
        &self,
        samples_names_to_run: HashSet<&str>,
        read_collection: &ReadCollection,
    ) -> Result<CoverageTable> {
        debug!("Sample names to run: {:?}", samples_names_to_run);
        let mut coverage_tables = Vec::with_capacity(4);

        let short_reads_to_map = read_collection.subset_short_reads(&samples_names_to_run);
        if short_reads_to_map.len() > 0 {
            info!("Mapping {} short reads.", short_reads_to_map.len());
            coverage_tables.push(self.run_coverm(short_reads_to_map, MappingMode::ShortRead)?);
        }

        let long_reads_to_map = read_collection.subset_long_reads(&samples_names_to_run);
        if long_reads_to_map.len() > 0 {
            info!("Mapping {} long reads.", long_reads_to_map.len());
            coverage_tables.push(self.run_coverm(long_reads_to_map, MappingMode::LongRead)?);
        }

        let short_bams_to_use = read_collection.subset_short_read_bams(&samples_names_to_run);
        if short_bams_to_use.len() > 0 {
            info!(
                "Calculating coverage for {} short read bams.",
                short_bams_to_use.len()
            );
            coverage_tables.push(self.run_coverm(short_bams_to_use, MappingMode::ShortBam)?);
        }

        let long_bams_to_use = read_collection.subset_long_read_bams(&samples_names_to_run);
        if long_bams_to_use.len() > 0 {
            info!(
                "Calculating coverage for {} long read bams.",
                long_bams_to_use.len()
            );
            coverage_tables.push(self.run_coverm(long_bams_to_use, MappingMode::LongBam)?);
        }

        CoverageTable::merge_many(coverage_tables)
    }

    fn run_coverm(
        &self,
        read_collection: ReadCollection,
        mode: MappingMode,
    ) -> Result<CoverageTable> {
        let mut coverm_command = Command::new("coverm");
        coverm_command
            .arg("contig")
            .arg("--threads")
            .arg(&format!("{}", self.threads))
            .arg("--min-covered-fraction")
            .arg(&format!("{}", self.trimming.min_covered_fraction));

        // Short and long reads go to separate CoverM invocations, so each needs its own
        // mapper. Unset means CoverM's own default, which is why the name is not validated
        // here: the list belongs to CoverM and drifts whenever CoverM adds one.
        let mapper = match mode {
            MappingMode::ShortRead => self.mapping.mapper.as_ref(),
            MappingMode::LongRead => self.mapping.longread_mapper.as_ref(),
            _ => None,
        };
        if matches!(mode, MappingMode::ShortRead | MappingMode::LongRead) {
            if let Some(mapper) = mapper {
                coverm_command.arg("--mapper").arg(mapper);
            }
            coverm_command.arg("--reference").arg(&self.assembly);
        }

        if let Some(minimap_params) = &self.mapping.minimap2_params {
            coverm_command.arg("--minimap2-params").arg(minimap_params);
        }
        if let Some(bwa_params) = &self.mapping.bwa_params {
            coverm_command.arg("--bwa-params").arg(bwa_params);
        }
        if let Some(min_read_aligned_length) = &self.filtering.min_read_aligned_length {
            coverm_command
                .arg("--min-read-aligned-length")
                .arg(&format!("{}", min_read_aligned_length));
        }
        if let Some(min_read_percent_identity) = &self.filtering.min_read_percent_identity {
            coverm_command
                .arg("--min-read-percent-identity")
                .arg(&format!("{}", min_read_percent_identity));
        }
        coverm_command
            .arg("--min-read-aligned-percent")
            .arg(&format!("{}", self.filtering.min_read_aligned_percent));
        if self.alignment.include_secondary {
            coverm_command.arg("--include-secondary");
        }
        if self.alignment.exclude_supplementary {
            coverm_command.arg("--exclude-supplementary");
        }

        coverm_command
            .arg("--contig-end-exclusion")
            .arg(&format!("{}", self.trimming.contig_end_exclusion))
            .arg("--trim-min")
            .arg(&format!("{}", self.trimming.trim_min))
            .arg("--trim-max")
            .arg(&format!("{}", self.trimming.trim_max));

        match mode {
            MappingMode::ShortRead | MappingMode::ShortBam => {
                coverm_command.arg("--methods").arg("metabat");
                if let Some(min_read_aligned_length_pair) =
                    &self.filtering.min_read_aligned_length_pair
                {
                    coverm_command
                        .arg("--min-read-aligned-length-pair")
                        .arg(&format!("{}", min_read_aligned_length_pair));
                }
                if let Some(min_read_percent_identity_pair) =
                    &self.filtering.min_read_percent_identity_pair
                {
                    coverm_command
                        .arg("--min-read-percent-identity-pair")
                        .arg(&format!("{}", min_read_percent_identity_pair));
                }
                if let Some(min_read_aligned_percent_pair) =
                    &self.filtering.min_read_aligned_percent_pair
                {
                    coverm_command
                        .arg("--min-read-aligned-percent-pair")
                        .arg(&format!("{}", min_read_aligned_percent_pair));
                }
                if self.alignment.proper_pairs_only {
                    coverm_command.arg("--proper-pairs-only");
                }
            }
            MappingMode::LongRead | MappingMode::LongBam => {
                coverm_command
                    .arg("--methods")
                    .arg("length")
                    .arg("trimmed_mean")
                    .arg("variance");
            }
        };

        read_collection.add_to_coverm_command(&mut coverm_command);

        let temp_file = tempfile::NamedTempFile::new()?;
        let temp_file_path = temp_file.path().to_str().unwrap();
        coverm_command.arg("--output-file").arg(temp_file_path);

        coverm_command
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped());

        match coverm_command.output() {
            Ok(output) => {
                if output.status.success() {
                    let coverage_table = CoverageTable::from_file(temp_file_path, mode)?;
                    Ok(coverage_table)
                } else {
                    Err(anyhow::anyhow!(
                        "Coverm failed with exit code: {} {} {}",
                        output.status,
                        std::str::from_utf8(output.stdout.as_slice())?,
                        std::str::from_utf8(output.stderr.as_slice())?
                    ))
                }
            }
            Err(e) => Err(anyhow::anyhow!("Coverm failed with error: {}", e)),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum MappingMode {
    ShortRead,
    LongRead,
    ShortBam,
    LongBam,
}

/// The version rosella's two table parsers were written against. Declared in `pixi.toml`
/// and `rosella.yml` as well, and asserted here so a stale environment fails at startup
/// rather than inside a parse.
const MINIMUM_COVERM: (u32, u32, u32) = (0, 6, 1);

fn check_coverm_is_installed() -> Result<()> {
    let reported = match Command::new("coverm").arg("--version").output() {
        Ok(output) if output.status.success() => {
            String::from_utf8_lossy(&output.stdout).into_owned()
        }
        Ok(output) => bail!(
            "`coverm --version` failed: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        ),
        Err(_) => bail!(
            "coverm is not on PATH. rosella maps reads and reads BAMs through it, so either \
             install it or pass a coverage table with --coverage-file"
        ),
    };

    let Some(version) = parse_version(&reported) else {
        debug!("Could not read a version out of `{}`", reported.trim());
        return Ok(());
    };
    if version < MINIMUM_COVERM {
        bail!(
            "coverm {}.{}.{} is older than the {}.{}.{} rosella's coverage table parsers \
             were written against",
            version.0,
            version.1,
            version.2,
            MINIMUM_COVERM.0,
            MINIMUM_COVERM.1,
            MINIMUM_COVERM.2
        );
    }
    debug!("coverm {}.{}.{}", version.0, version.1, version.2);
    Ok(())
}

fn parse_version(reported: &str) -> Option<(u32, u32, u32)> {
    let mut parts = reported.split_whitespace().nth(1)?.split('.');
    Some((
        parts.next()?.parse().ok()?,
        parts.next()?.parse().ok()?,
        parts.next().unwrap_or("0").parse().unwrap_or(0),
    ))
}
