use std::{collections::HashSet, process::Command};

use anyhow::Result;
use log::{info, debug};

use crate::coverage::{coverage_calculator::ReadCollection, coverage_table::CoverageTable};


pub struct CovermEngine<'a> {
    assembly: &'a str,
    threads: usize,
    settings: CovermSettings,
}

impl<'a> CovermEngine<'a> {
    pub fn new(m: &'a clap::ArgMatches) -> Result<Self> {
        // get the assembly path
        let assembly = m.get_one::<String>("assembly").unwrap();

        // create output directory
        let output_directory = m.get_one::<String>("output-directory").unwrap();
        // create the output directory but do not fail if it already exists
        std::fs::create_dir_all(&output_directory)?;

        let settings = CovermSettings::new(m)?;

        Ok(
            Self {
                assembly,
                threads: *m.get_one::<usize>("threads").unwrap(),
                settings
            }
        )
    }

    pub fn run(&self, samples_names_to_run: HashSet<&str>, read_collection: &ReadCollection) -> Result<CoverageTable> {
        
        debug!("Sample names to run: {:?}", samples_names_to_run);
        // collect short reads that need to be mapped
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
            info!("Calculating coverage for {} short read bams.", short_bams_to_use.len());
            coverage_tables.push(self.run_coverm(short_bams_to_use, MappingMode::ShortBam)?);
        }

        let long_bams_to_use = read_collection.subset_long_read_bams(&samples_names_to_run);
        if long_bams_to_use.len() > 0 {
            info!("Calculating coverage for {} long read bams.", long_bams_to_use.len());
            coverage_tables.push(self.run_coverm(long_bams_to_use, MappingMode::LongBam)?);
        }

        CoverageTable::merge_many(coverage_tables)
    }

    fn run_coverm(&self, read_collection: ReadCollection, mode: MappingMode) -> Result<CoverageTable> {
        let mut coverm_command = Command::new("coverm");
        coverm_command
            .arg("contig")
            .arg("--threads").arg(&format!("{}", self.threads))
            .arg("--min-covered-fraction").arg(&format!("{}", self.settings.min_covered_fraction));

        // add the mapper and reference
        match mode {
            MappingMode::ShortRead => {
                // map short reads
                coverm_command
                    .arg("--mapper")
                    .arg(&self.settings.short_read_mapper)
                    .arg("--reference")
                    .arg(&self.assembly);
            }
            MappingMode::LongRead => {
                // map long reads
                coverm_command
                    .arg("--mapper")
                    .arg(&self.settings.long_read_mapper)
                    .arg("--reference")
                    .arg(&self.assembly);
            }
            _ => {
                // nothing yet
            }
        };

        if let Some(minimap_params) = &self.settings.minimap_params {
            coverm_command
                .arg("--minimap2-params")
                .arg(minimap_params);
        }
        if let Some(bwa_params) = &self.settings.bwa_params {
            coverm_command
                .arg("--bwa-params")
                .arg(bwa_params);
        }
        if let Some(min_read_aligned_length) = &self.settings.min_read_aligned_length {
            coverm_command
                .arg("--min-read-aligned-length")
                .arg(&format!("{}", min_read_aligned_length));
        }
        if let Some(min_read_percent_identity) = &self.settings.min_read_percent_identity {
            coverm_command
                .arg("--min-read-percent-identity")
                .arg(&format!("{}", min_read_percent_identity));
        }
        if let Some(min_read_aligned_percent) = &self.settings.min_read_aligned_percent {
            coverm_command
                .arg("--min-read-aligned-percent")
                .arg(&format!("{}", min_read_aligned_percent));
        }
        if self.settings.include_secondary {
            coverm_command
                .arg("--include-secondary");
        }
        if self.settings.exclude_supplementary {
            coverm_command
                .arg("--exclude-supplementary");
        }

        coverm_command
            .arg("--contig-end-exclusion")
            .arg(&format!("{}", self.settings.contig_end_exclusion))
            .arg("--trim-min")
            .arg(&format!("{}", self.settings.trim_min))
            .arg("--trim-max")
            .arg(&format!("{}", self.settings.trim_max));

        // add the short or long read read filter settings.
        match mode {
            MappingMode::ShortRead | MappingMode::ShortBam => {
                // map short reads
                coverm_command
                    .arg("--methods")
                    .arg("metabat");
                if let Some(min_read_aligned_length_pair) = &self.settings.min_read_aligned_length_pair {
                    coverm_command
                        .arg("--min-read-aligned-length-pair")
                        .arg(&format!("{}", min_read_aligned_length_pair));
                }
                if let Some(min_read_percent_identity_pair) = &self.settings.min_read_percent_identity_pair {
                    coverm_command
                        .arg("--min-read-percent-identity-pair")
                        .arg(&format!("{}", min_read_percent_identity_pair));
                }
                if let Some(min_read_aligned_percent_pair) = &self.settings.min_read_aligned_percent_pair {
                    coverm_command
                        .arg("--min-read-aligned-percent-pair")
                        .arg(&format!("{}", min_read_aligned_percent_pair));
                }
                if self.settings.proper_pairs_only {
                    coverm_command
                        .arg("--proper-pairs-only");
                }

            }
            MappingMode::LongRead | MappingMode::LongBam => {
                // need to be less stringent with read matching
                coverm_command
                    .arg("--methods")
                    .arg("length")
                    .arg("trimmed_mean")
                    .arg("variance");
                
            }
        };

        // add in the read collection
        read_collection.add_to_coverm_command(&mut coverm_command);

        // write to temp file
        let temp_file = tempfile::NamedTempFile::new()?;
        let temp_file_path = temp_file.path().to_str().unwrap();
        coverm_command
            .arg("--output-file")
            .arg(temp_file_path);

        // pipe the stdout and stderr
        coverm_command
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped());

        match coverm_command.output() {
            Ok(output) => {
                if output.status.success() {
                    // parse the output
                    let coverage_table = CoverageTable::from_file(temp_file_path, mode)?;
                    Ok(coverage_table)
                } else {
                    Err(anyhow::anyhow!("Coverm failed with exit code: {} {} {}", 
                        output.status, 
                        std::str::from_utf8(output.stdout.as_slice())?, 
                        std::str::from_utf8(output.stderr.as_slice())?
                    ))
                }
            }
            Err(e) => {
                Err(anyhow::anyhow!("Coverm failed with error: {}", e))
            }
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

struct CovermSettings {
    short_read_mapper: String,
    long_read_mapper: String,
    minimap_params: Option<String>,
    bwa_params: Option<String>,
    min_read_aligned_length: Option<u32>,
    min_read_percent_identity: Option<f32>,
    min_read_aligned_percent: Option<f32>,
    min_read_aligned_length_pair: Option<u32>,
    min_read_percent_identity_pair: Option<f32>,
    min_read_aligned_percent_pair: Option<f32>,
    proper_pairs_only: bool,
    min_covered_fraction: f32,
    include_secondary: bool,
    exclude_supplementary: bool,
    contig_end_exclusion: usize,
    trim_min: f32,
    trim_max: f32,
}

impl CovermSettings {
    pub fn new(m: &clap::ArgMatches) -> Result<Self> {
        let minimap_params = m.get_one::<String>("minimap2-params").cloned();
        let bwa_params = m.get_one::<String>("bwa-params").cloned();
        let min_read_aligned_length = m.get_one::<u32>("min-read-aligned-length").cloned();
        let min_read_percent_identity = m.get_one::<f32>("min-read-percent-identity").cloned();
        let min_read_aligned_percent = m.get_one::<f32>("min-read-aligned-percent").cloned();
        let min_read_aligned_length_pair = m.get_one::<u32>("min-read-aligned-length-pair").cloned();
        let min_read_percent_identity_pair = m.get_one::<f32>("min-read-percent-identity-pair").cloned();
        let min_read_aligned_percent_pair = m.get_one::<f32>("min-read-aligned-percent-pair").cloned();
        let proper_pairs_only = m.get_flag("proper-pairs-only");
        let min_covered_fraction = m.get_one::<f32>("min-covered-fraction").unwrap().clone();
        let include_secondary = m.get_flag("include-secondary");
        let exclude_supplementary = m.get_flag("exclude-supplementary");
        let contig_end_exclusion = m.get_one::<usize>("contig-end-exclusion").unwrap().clone();
        let trim_min = m.get_one::<f32>("trim-min").unwrap().clone();
        let trim_max = m.get_one::<f32>("trim-max").unwrap().clone();

        Ok(
            Self {
                short_read_mapper: m.get_one::<String>("mapper").unwrap().clone(),
                long_read_mapper: m.get_one::<String>("longread-mapper").unwrap().clone(),
                minimap_params,
                bwa_params,
                min_read_aligned_length,
                min_read_percent_identity,
                min_read_aligned_percent,
                min_read_aligned_length_pair,
                min_read_percent_identity_pair,
                min_read_aligned_percent_pair,
                proper_pairs_only,
                min_covered_fraction,
                include_secondary,
                exclude_supplementary,
                contig_end_exclusion,
                trim_min,
                trim_max,
            }
        )
    }
}