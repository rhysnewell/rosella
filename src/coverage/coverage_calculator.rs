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
    /// Only the mapping paths need it, and `refine` can be given both tables instead.
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
        // find previously calculated samples
        let previous_sample_names = self.find_previous_calculated_samples()?;
        debug!("previous sample names: {:?}", previous_sample_names);

        match (previous_sample_names, &self.read_collection) {
            (Some(previous_samples), Some(read_collection)) => {
                // if there are previously calculated samples, then we want to check if the
                // samples we are calculating now are already present in the previous samples
                // if they are, then we want to skip them

                let mut samples_to_calculate = HashSet::new();
                for sample_name in read_collection.sample_names() {
                    if !previous_samples.contains(sample_name) {
                        samples_to_calculate.insert(sample_name);
                    }
                }

                if samples_to_calculate.is_empty() {
                    // if there are no samples to calculate, then we want to return the
                    // previously calculated coverage table
                    let coverage_table =
                        CoverageTable::from_any_file(self.coverage_table_path.as_ref().unwrap())?;
                    return Ok(coverage_table);
                }
                let coverm_engine = CovermEngine::new(inputs)?;
                let new_coverages = coverm_engine.run(samples_to_calculate, read_collection)?;

                match &self.coverage_table_path {
                    Some(old) => {
                        // merge old and new coverages
                        let mut old_coverages = CoverageTable::from_any_file(old)?;
                        old_coverages.merge(new_coverages)?;
                        let output_file = format!("{}/coverage.tsv", self.output_directory);
                        old_coverages.write(output_file)?;
                        Ok(old_coverages)
                    }
                    None => Ok(new_coverages),
                }
            }
            (None, Some(read_collection)) => {
                // if there are no previously calculated samples, then we want to run coverm
                // on all samples
                let sample_names = read_collection
                    .sample_names()
                    .into_iter()
                    .collect::<HashSet<_>>();
                let coverm_engine = CovermEngine::new(inputs)?;
                let mut coverages = coverm_engine.run(sample_names, read_collection)?;
                let output_file = format!("{}/coverage.tsv", self.output_directory);
                coverages.write(output_file)?;
                Ok(coverages)
            }
            (Some(_), None) => {
                // if there are previously calculated samples, but no reads, then we want to
                // return the previously calculated coverage table
                let coverage_table =
                    CoverageTable::from_any_file(self.coverage_table_path.as_ref().unwrap())?;
                Ok(coverage_table)
            }
            (None, None) => Err(anyhow!("No coverage file or reads provided.")),
        }
    }

    /// If a coverage table is present, we want to check what samples are present
    /// in the header. This function will return a set of sample names.
    fn find_previous_calculated_samples(&self) -> Result<Option<HashSet<String>>> {
        match &self.coverage_table_path {
            Some(path) => Ok(Some(
                CoverageTable::sample_names_in(path)?.into_iter().collect(),
            )),
            None => Ok(None),
        }
    }
}

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
        let mut interleaved: Option<Vec<_>> = None;
        let mut unpaired: Option<Vec<_>> = None;

        if let Some(read1_paths) = &self.forward_read_paths {
            let mut inner_read1 = Vec::new();
            let mut inner_read2 = Vec::new();
            for (read1_path, read2_path) in read1_paths
                .iter()
                .zip(self.reverse_read_paths.as_ref().unwrap())
            {
                let read1_path = Path::new(read1_path);
                let read2_path = Path::new(read2_path);
                let sample_name = read1_path.file_name().unwrap().to_str().unwrap();
                if sample_names_to_map.contains(sample_name) {
                    inner_read1.push(read1_path.to_str().unwrap().to_string());
                    inner_read2.push(read2_path.to_str().unwrap().to_string());
                }
            }
            read1 = non_empty(inner_read1);
            read2 = non_empty(inner_read2);
        }

        if let Some(interleaved_paths) = &self.interleaved_read_paths {
            let mut inner_interleaved = Vec::new();
            for interleaved_path in interleaved_paths {
                let interleaved_path = Path::new(interleaved_path);
                let sample_name = interleaved_path.file_name().unwrap().to_str().unwrap();
                if sample_names_to_map.contains(sample_name) {
                    inner_interleaved.push(interleaved_path.to_str().unwrap().to_string());
                }
            }
            interleaved = non_empty(inner_interleaved);
        }

        if let Some(unpaired_paths) = &self.unpaired_read_paths {
            let mut inner_unpaired = Vec::new();
            for unpaired_path in unpaired_paths {
                let unpaired_path = Path::new(unpaired_path);
                let sample_name = unpaired_path.file_name().unwrap().to_str().unwrap();
                if sample_names_to_map.contains(sample_name) {
                    inner_unpaired.push(unpaired_path.to_str().unwrap().to_string());
                }
            }
            unpaired = non_empty(inner_unpaired);
        }

        Self {
            forward_read_paths: read1,
            reverse_read_paths: read2,
            interleaved_read_paths: interleaved,
            unpaired_read_paths: unpaired,
            long_read_paths: None,
            short_read_bam_paths: None,
            long_read_bam_paths: None,
        }
    }

    /// returns a ReadCollection that only contains the
    /// sample names specified in the sample_names_to_map
    /// subset to only longread samples and no BAM files
    pub fn subset_long_reads(&self, sample_names_to_map: &HashSet<&str>) -> Self {
        let mut long_reads: Option<Vec<_>> = None;

        if let Some(long_read_paths) = &self.long_read_paths {
            debug!("long read paths: {:?}", long_read_paths);
            let mut inner_long_reads = Vec::new();
            for long_read_path in long_read_paths {
                let long_read_path = Path::new(long_read_path);
                let sample_name = long_read_path.file_name().unwrap().to_str().unwrap();
                if sample_names_to_map.contains(sample_name) {
                    inner_long_reads.push(long_read_path.to_str().unwrap().to_string());
                }
            }
            long_reads = non_empty(inner_long_reads);
        }

        Self {
            forward_read_paths: None,
            reverse_read_paths: None,
            interleaved_read_paths: None,
            unpaired_read_paths: None,
            long_read_paths: long_reads,
            short_read_bam_paths: None,
            long_read_bam_paths: None,
        }
    }

    /// returns a ReadCollection that only contains the
    /// sample names specified in the sample_names_to_map
    /// subset to only shortread BAM files
    pub fn subset_short_read_bams(&self, sample_names_to_keep: &HashSet<&str>) -> Self {
        let mut short_read_bams: Option<Vec<_>> = None;

        if let Some(short_read_bam_paths) = &self.short_read_bam_paths {
            let mut inner_short_read_bams = Vec::new();
            for short_read_bam_path in short_read_bam_paths {
                let short_read_bam_path = Path::new(short_read_bam_path);
                let sample_name = short_read_bam_path.file_name().unwrap().to_str().unwrap();
                if sample_names_to_keep.contains(sample_name) {
                    inner_short_read_bams.push(short_read_bam_path.to_str().unwrap().to_string());
                }
            }
            short_read_bams = non_empty(inner_short_read_bams);
        }

        Self {
            forward_read_paths: None,
            reverse_read_paths: None,
            interleaved_read_paths: None,
            unpaired_read_paths: None,
            long_read_paths: None,
            short_read_bam_paths: short_read_bams,
            long_read_bam_paths: None,
        }
    }

    /// returns a ReadCollection that only contains the
    /// sample names specified in the sample_names_to_map
    /// subset to only longread BAM files
    pub fn subset_long_read_bams(&self, sample_names_to_keep: &HashSet<&str>) -> Self {
        let mut long_read_bams: Option<Vec<_>> = None;

        if let Some(long_read_bam_paths) = &self.long_read_bam_paths {
            let mut inner_long_read_bams = Vec::new();
            for long_read_bam_path in long_read_bam_paths {
                let long_read_bam_path = Path::new(long_read_bam_path);
                let sample_name = long_read_bam_path.file_name().unwrap().to_str().unwrap();
                if sample_names_to_keep.contains(sample_name) {
                    inner_long_read_bams.push(long_read_bam_path.to_str().unwrap().to_string());
                }
            }
            long_read_bams = non_empty(inner_long_read_bams);
        }

        Self {
            forward_read_paths: None,
            reverse_read_paths: None,
            interleaved_read_paths: None,
            unpaired_read_paths: None,
            long_read_paths: None,
            short_read_bam_paths: None,
            long_read_bam_paths: long_read_bams,
        }
    }

    /// returns a vector of all of the sample names
    /// in the read collection.
    /// The sample names are the file stem of the read file
    /// So no path information is included in the sample name
    pub fn sample_names(&self) -> Vec<&str> {
        let mut sample_names = Vec::new();

        // for read1 and read2, we only need read1 names
        if let Some(read1) = &self.forward_read_paths {
            for path in read1 {
                let path = Path::new(path);
                sample_names.push(path.file_name().unwrap().to_str().unwrap());
            }
        }

        if let Some(interleaved) = &self.interleaved_read_paths {
            for path in interleaved {
                let path = Path::new(path);
                sample_names.push(path.file_name().unwrap().to_str().unwrap());
            }
        }

        if let Some(unpaired) = &self.unpaired_read_paths {
            for path in unpaired {
                let path = Path::new(path);
                sample_names.push(path.file_name().unwrap().to_str().unwrap());
            }
        }

        if let Some(long_reads) = &self.long_read_paths {
            for path in long_reads {
                let path = Path::new(path);
                sample_names.push(path.file_name().unwrap().to_str().unwrap());
            }
        }

        if let Some(short_read_bams) = &self.short_read_bam_paths {
            for path in short_read_bams {
                let path = Path::new(path);
                sample_names.push(path.file_name().unwrap().to_str().unwrap());
            }
        }

        if let Some(long_read_bams) = &self.long_read_bam_paths {
            for path in long_read_bams {
                let path = Path::new(path);
                sample_names.push(path.file_name().unwrap().to_str().unwrap());
            }
        }

        debug!("sample names: {:?}", sample_names);

        sample_names
    }

    pub fn len(&self) -> usize {
        let mut len = 0;

        if let Some(read1) = &self.forward_read_paths {
            len += read1.len();
        }

        if let Some(interleaved) = &self.interleaved_read_paths {
            len += interleaved.len();
        }

        if let Some(unpaired) = &self.unpaired_read_paths {
            len += unpaired.len();
        }

        if let Some(long_reads) = &self.long_read_paths {
            len += long_reads.len();
        }

        if let Some(short_read_bams) = &self.short_read_bam_paths {
            len += short_read_bams.len();
        }

        if let Some(long_read_bams) = &self.long_read_bam_paths {
            len += long_read_bams.len();
        }

        len
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
fn non_empty(paths: Vec<String>) -> Option<Vec<String>> {
    (!paths.is_empty()).then_some(paths)
}
