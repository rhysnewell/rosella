use std::{collections::HashSet, process::Command};

use anyhow::Result;
use bird_tool_utils::command::finish_command_safely;
use log::{warn, info};

use crate::coverage::{coverage_calculator::ReadCollection, coverage_table::CoverageTable};


pub struct CovermEngine<'a> {
    assembly: &'a str,
    output_directory: &'a str,
    threads: usize,
    settings: CovermSettings,
}

impl<'a> CovermEngine<'a> {
    pub fn new(m: &clap::ArgMatches) -> Result<Self> {
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
                output_directory,
                threads: *m.get_one::<usize>("threads").unwrap(),
                settings
            }
        )
    }

    pub fn run(&self, samples_names_to_run: HashSet<&str>, read_collection: &Option<ReadCollection>) -> Result<Option<CoverageTable>> {
        // collect short reads that need to be mapped
        let coverage_table = match read_collection {
            Some(read_collection) => {
                let mut coverage_tables = Vec::with_capacity(4);

                let short_reads_to_map = read_collection.subset_short_reads(&samples_names_to_run);
                info!("Mapping {} short reads.", short_reads_to_map.len());
                coverage_tables.push(self.run_coverm(short_reads_to_map, MappingMode::ShortRead)?);

                let long_reads_to_map = read_collection.subset_long_reads(&samples_names_to_run);
                info!("Mapping {} long reads.", long_reads_to_map.len());
                coverage_tables.push(self.run_coverm(long_reads_to_map, MappingMode::LongRead)?);

                let short_bams_to_use = read_collection.subset_short_read_bams(&samples_names_to_run);
                info!("Calculating coverage for {} short read bams.", short_bams_to_use.len());
                coverage_tables.push(self.run_coverm(short_bams_to_use, MappingMode::ShortBam)?);

                let long_bams_to_use = read_collection.subset_long_read_bams(&samples_names_to_run);
                info!("Calculating coverage for {} long read bams.", long_bams_to_use.len());
                coverage_tables.push(self.run_coverm(long_bams_to_use, MappingMode::LongBam)?);

                Some(CoverageTable::merge(coverage_tables))
            }
            None => {
                warn!("No ReadCollection provided. Skipping mapping and coverage calculation.");
                None
            }
        };


        Ok(coverage_table)
    }

    fn run_coverm(&self, read_collection: ReadCollection, mode: MappingMode) -> Result<CoverageTable> {
        let mut coverm_command = Command::new("coverm");
        coverm_command
            .arg("contig")
            .arg("--threads").arg(&format!("{}", self.threads));

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
                    .arg(&self.assembly)
                    .arg("--threads");
            }
            _ => {
                // nothing yet
            }
        };

        // add the short or long read read filter settings.
        match mode {
            MappingMode::ShortRead | MappingMode::ShortBam => {
                // map short reads
                coverm_command
                    .arg("--mode")
                    .arg("metabat");
            }
            MappingMode::LongRead | MappingMode::LongBam => {
                // need to be less stringent with read matching
                coverm_command
                    .arg("--mode")
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
            .arg("--output")
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
                    Err(anyhow::anyhow!("Coverm failed with exit code: {}", output.status))
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
}

impl CovermSettings {
    pub fn new(m: &clap::ArgMatches) -> Result<Self> {
        Ok(
            Self {
                short_read_mapper: m.get_one::<String>("mapper").unwrap().clone(),
                long_read_mapper: m.get_one::<String>("longread-mapper").unwrap().clone(),
            }
        )
    }
}