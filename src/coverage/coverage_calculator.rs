use anyhow::{Result, anyhow};
use log::debug;
use std::{collections::HashSet, path::Path, process::Command};

use super::coverage_table::CoverageTable;
use crate::cli::{AlignmentFlags, CoverageSource, CoverageTrimming, MappingParams, ReadFiltering};
use crate::external::coverm_engine::CovermEngine;

/// Everything a coverage table needs, gathered from whichever subcommand asked for one.
pub struct CoverageInputs<'a> {
    pub assembly: Option<&'a str>,
    pub output_directory: &'a str,
    pub threads: usize,
    pub coverage: &'a CoverageSource,
    pub mapping: &'a MappingParams,
    pub filtering: &'a ReadFiltering,
    pub alignment: &'a AlignmentFlags,
    pub trimming: &'a CoverageTrimming,
}

impl<'a> CoverageInputs<'a> {
    /// Only the mapping paths need it, and `score` never maps.
    pub fn assembly(&self) -> Result<&'a str> {
        self.assembly
            .ok_or_else(|| anyhow!("mapping reads needs an assembly, so pass --assembly"))
    }
}

/// Coverage is either calculated from the reads through CoverM or read from a table.
pub fn calculate_coverage(inputs: &CoverageInputs) -> Result<CoverageTable> {
    let mut engine = CoverageCalculatorEngine::new(inputs)?;
    engine.run(inputs)
}

struct CoverageCalculatorEngine {
    read_collection: Option<ReadCollection>,
    coverage_table_path: Option<String>,
    output_directory: String,
}

impl CoverageCalculatorEngine {
    pub fn new(inputs: &CoverageInputs) -> Result<Self> {
        let output_directory = inputs.output_directory.to_string();
        std::fs::create_dir_all(&output_directory)?;

        let coverage_table_path = match &inputs.coverage.coverage_file {
            Some(coverage_table_path) => Some(coverage_table_path.clone()),
            None => {
                let coverage_table_path = format!("{}/coverage.tsv", output_directory);
                if std::path::Path::new(&coverage_table_path).exists() {
                    Some(coverage_table_path)
                } else {
                    None
                }
            }
        };

        let read_collection = Some(ReadCollection::new(inputs.coverage)?);

        Ok(Self {
            read_collection,
            coverage_table_path,
            output_directory,
        })
    }

    pub fn run(&mut self, inputs: &CoverageInputs) -> Result<CoverageTable> {
        let previous_sample_names = self.find_previous_calculated_samples()?;
        debug!("previous sample names: {:?}", previous_sample_names);

        match (previous_sample_names, &self.read_collection) {
            (Some(previous_samples), Some(read_collection)) => {
                let mut samples_to_calculate = HashSet::new();
                for sample_name in read_collection.sample_names() {
                    if !previous_samples.contains(sample_name) {
                        samples_to_calculate.insert(sample_name);
                    }
                }

                if samples_to_calculate.is_empty() {
                    let coverage_table =
                        CoverageTable::from_any_file(self.coverage_table_path.as_ref().unwrap())?;
                    return Ok(coverage_table);
                }
                let coverm_engine = CovermEngine::new(inputs)?;
                let new_coverages = coverm_engine.run(samples_to_calculate, read_collection)?;

                match &self.coverage_table_path {
                    Some(old) => {
                        let mut old_coverages = CoverageTable::from_any_file(old)?;
                        old_coverages.merge(new_coverages)?;
                        old_coverages.align_to(&read_collection.sample_names());
                        let output_file = format!("{}/coverage.tsv", self.output_directory);
                        old_coverages.write(output_file)?;
                        Ok(old_coverages)
                    }
                    None => Ok(new_coverages),
                }
            }
            (None, Some(read_collection)) => {
                let sample_names = read_collection
                    .sample_names()
                    .into_iter()
                    .collect::<HashSet<_>>();
                let coverm_engine = CovermEngine::new(inputs)?;
                let coverages = coverm_engine.run(sample_names, read_collection)?;
                let output_file = format!("{}/coverage.tsv", self.output_directory);
                coverages.write(output_file)?;
                Ok(coverages)
            }
            (Some(_), None) => {
                let coverage_table =
                    CoverageTable::from_any_file(self.coverage_table_path.as_ref().unwrap())?;
                Ok(coverage_table)
            }
            (None, None) => Err(anyhow!("No coverage file or reads provided.")),
        }
    }

    fn find_previous_calculated_samples(&self) -> Result<Option<HashSet<String>>> {
        match &self.coverage_table_path {
            Some(path) => Ok(Some(
                CoverageTable::sample_names_in(path)?.into_iter().collect(),
            )),
            None => Ok(None),
        }
    }
}

#[derive(Default)]
pub struct ReadCollection {
    forward_read_paths: Option<Vec<String>>,
    reverse_read_paths: Option<Vec<String>>,
    interleaved_read_paths: Option<Vec<String>>,
    unpaired_read_paths: Option<Vec<String>>,
    long_read_paths: Option<Vec<String>>,
    short_read_bam_paths: Option<Vec<String>>,
    long_read_bam_paths: Option<Vec<String>>,
}

impl ReadCollection {
    pub fn new(source: &CoverageSource) -> Result<Self> {
        if source.read1.len() != source.read2.len() {
            bail!(
                "When specifying paired reads with the -1 and -2 flags, there must be equal \
                 numbers specified. Instead found {} and {} respectively",
                source.read1.len(),
                source.read2.len()
            );
        }
        if source.coupled.len() % 2 != 0 {
            bail!(
                "The --coupled flag must be set with pairs of read sets, but an odd number \
                 ({}) was specified",
                source.coupled.len()
            );
        }

        let mut read1 = source.read1.clone();
        let mut read2 = source.read2.clone();
        for pair in source.coupled.chunks(2) {
            read1.push(pair[0].clone());
            read2.push(pair[1].clone());
        }

        Ok(Self {
            forward_read_paths: non_empty(read1),
            reverse_read_paths: non_empty(read2),
            interleaved_read_paths: non_empty(source.interleaved.clone()),
            unpaired_read_paths: non_empty(source.single.clone()),
            long_read_paths: non_empty(source.longreads.clone()),
            short_read_bam_paths: non_empty(source.bam_files.clone()),
            long_read_bam_paths: non_empty(source.longread_bam_files.clone()),
        })
    }

    pub fn subset_short_reads(&self, sample_names_to_map: &HashSet<&str>) -> Self {
        let mut read1: Option<Vec<_>> = None;
        let mut read2: Option<Vec<_>> = None;

        if let Some(read1_paths) = &self.forward_read_paths {
            let mut inner_read1 = Vec::new();
            let mut inner_read2 = Vec::new();
            for (read1_path, read2_path) in read1_paths
                .iter()
                .zip(self.reverse_read_paths.as_ref().unwrap())
            {
                if sample_names_to_map.contains(sample_of(read1_path)) {
                    inner_read1.push(read1_path.clone());
                    inner_read2.push(read2_path.clone());
                }
            }
            read1 = non_empty(inner_read1);
            read2 = non_empty(inner_read2);
        }

        Self {
            forward_read_paths: read1,
            reverse_read_paths: read2,
            interleaved_read_paths: kept(&self.interleaved_read_paths, sample_names_to_map),
            unpaired_read_paths: kept(&self.unpaired_read_paths, sample_names_to_map),
            ..Self::default()
        }
    }

    pub fn subset_long_reads(&self, sample_names_to_map: &HashSet<&str>) -> Self {
        Self {
            long_read_paths: kept(&self.long_read_paths, sample_names_to_map),
            ..Self::default()
        }
    }

    pub fn subset_short_read_bams(&self, sample_names_to_keep: &HashSet<&str>) -> Self {
        Self {
            short_read_bam_paths: kept(&self.short_read_bam_paths, sample_names_to_keep),
            ..Self::default()
        }
    }

    pub fn subset_long_read_bams(&self, sample_names_to_keep: &HashSet<&str>) -> Self {
        Self {
            long_read_bam_paths: kept(&self.long_read_bam_paths, sample_names_to_keep),
            ..Self::default()
        }
    }

    fn arms(&self) -> [&Option<Vec<String>>; 6] {
        [
            &self.forward_read_paths,
            &self.interleaved_read_paths,
            &self.unpaired_read_paths,
            &self.long_read_paths,
            &self.short_read_bam_paths,
            &self.long_read_bam_paths,
        ]
    }

    /// Reverse reads are left out: a pair is one sample, named for its forward file.
    pub fn sample_names(&self) -> Vec<&str> {
        let sample_names = self
            .arms()
            .into_iter()
            .flatten()
            .flatten()
            .map(|path| sample_of(path))
            .collect::<Vec<_>>();
        debug!("sample names: {:?}", sample_names);
        sample_names
    }

    pub fn len(&self) -> usize {
        self.arms()
            .into_iter()
            .map(|arm| arm.as_ref().map_or(0, Vec::len))
            .sum()
    }

    pub fn add_to_coverm_command(&self, coverm_command: &mut Command) {
        if let Some(read1) = &self.forward_read_paths {
            coverm_command.arg("-1");
            for path in read1 {
                coverm_command.arg(path);
            }
        }

        if let Some(read2) = &self.reverse_read_paths {
            coverm_command.arg("-2");
            for path in read2 {
                coverm_command.arg(path);
            }
        }

        if let Some(interleaved) = &self.interleaved_read_paths {
            coverm_command.arg("--interleaved");
            for path in interleaved {
                coverm_command.arg(path);
            }
        }

        if let Some(unpaired) = &self.unpaired_read_paths {
            coverm_command.arg("--single");
            for path in unpaired {
                coverm_command.arg(path);
            }
            return;
        }

        if let Some(long_reads) = &self.long_read_paths {
            coverm_command.arg("--single");
            for path in long_reads {
                coverm_command.arg(path);
            }
            return;
        }

        if let Some(short_read_bams) = &self.short_read_bam_paths {
            coverm_command.arg("-b");
            for path in short_read_bams {
                coverm_command.arg(path);
            }
            return;
        }

        if let Some(long_read_bams) = &self.long_read_bam_paths {
            coverm_command.arg("-b");
            for path in long_read_bams {
                coverm_command.arg(path);
            }
        }
    }
}

/// An empty list still emitted a bare `-1` or `--single` with nothing after it, which
/// CoverM reads as the next flag's value.
fn sample_of(path: &str) -> &str {
    Path::new(path).file_name().unwrap().to_str().unwrap()
}

fn kept(paths: &Option<Vec<String>>, wanted: &HashSet<&str>) -> Option<Vec<String>> {
    non_empty(
        paths
            .iter()
            .flatten()
            .filter(|path| wanted.contains(sample_of(path)))
            .cloned()
            .collect(),
    )
}

fn non_empty(paths: Vec<String>) -> Option<Vec<String>> {
    (!paths.is_empty()).then_some(paths)
}
