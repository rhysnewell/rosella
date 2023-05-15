use std::{collections::HashSet, path::Path, process::Command};
use anyhow::{Result, anyhow};
use hnsw_rs::prelude::Distance;
use itertools::izip;
use ndarray::{ArrayView, prelude::*};
use statrs::distribution::{Normal, ContinuousCDF};

use crate::external::coverm_engine::{CovermEngine, MappingMode};
use super::coverage_table::CoverageTable;


/// This module contains the coverage calculator logic.
/// Coverage are either calculate from the reads using CoverM
/// or retrieved from a pre-calculated file.
pub fn calculate_coverage(m: &clap::ArgMatches) -> Result<CoverageTable> {
    let mut coverm_engine = CoverageCalculatorEngine::new(m)?;
    let coverage_table = coverm_engine.run(m)?;
    Ok(coverage_table)
}

struct CoverageCalculatorEngine {
    read_collection: Option<ReadCollection>,
    coverage_table_path: Option<String>,
    output_directory: String,
}

impl CoverageCalculatorEngine {
    pub fn new(m: &clap::ArgMatches) -> Result<Self> {
        // create output directory
        let output_directory = m.get_one::<String>("output-directory").unwrap().clone();
        // create the output directory but do not fail if it already exists
        std::fs::create_dir_all(&output_directory)?;

        // check if coverage file is provided or exists in output directory
        let coverage_table_path = match m.get_one::<String>("coverage-file") {
            Some(coverage_table_path) => {
                Some(coverage_table_path.clone())
            }
            None => {
                let coverage_table_path = format!("{}/coverage.tsv", output_directory);
                if std::path::Path::new(&coverage_table_path).exists() {
                    Some(coverage_table_path)
                } else {
                    None
                }
            }
        };

        // if coverage table is none, then check for the reads
        let read_collection = match coverage_table_path {
            Some(_) => None,
            None => {
                Some(ReadCollection::new(m)?)
            }
        };

        Ok(
            Self {
                read_collection,
                coverage_table_path,
                output_directory,
            }
        )
    }

    pub fn run(&mut self, m: &clap::ArgMatches) -> Result<CoverageTable> {
        // find previously calculated samples
        let previous_sample_names = self.find_previous_calculated_samples()?;
        
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

                let coverm_engine = CovermEngine::new(m)?;
                let new_coverages = coverm_engine.run(samples_to_calculate, read_collection)?;

                match &self.coverage_table_path {
                    Some(old) => {
                        // merge old and new coverages
                        let mut old_coverages = CoverageTable::from_file(&old, MappingMode::ShortBam)?;
                        old_coverages.merge(new_coverages)?;
                        let output_file = format!("{}/coverage.tsv", self.output_directory);
                        old_coverages.write(output_file)?;
                        Ok(old_coverages)
                    },
                    None => {
                        Ok(new_coverages)
                    }
                }
            },
            (None, Some(read_collection)) => {
                // if there are no previously calculated samples, then we want to run coverm
                // on all samples
                let sample_names = read_collection.sample_names().into_iter().collect::<HashSet<_>>();
                let coverm_engine = CovermEngine::new(m)?;
                let coverages = coverm_engine.run(sample_names, read_collection)?;
                let output_file = format!("{}/coverage.tsv", self.output_directory);
                coverages.write(output_file)?;
                Ok(coverages)
            },
            (Some(_), None) => {
                // if there are previously calculated samples, but no reads, then we want to
                // return the previously calculated coverage table
                let coverage_table = CoverageTable::from_file(&self.coverage_table_path.as_ref().unwrap(), MappingMode::ShortBam)?;
                Ok(coverage_table)
            },
            (None, None) => {
                Err(anyhow!("No coverage file or reads provided."))
            }
        }
    }

    /// If a coverage table is present, we want to check what samples are present
    /// in the header. This function will return a set of sample names.
    fn find_previous_calculated_samples(&self) -> Result<Option<HashSet<String>>> {
        match &self.coverage_table_path {
            Some(coverage_table_path) => {
                let mut previous_sample_names = HashSet::new();
                let mut reader = csv::ReaderBuilder::new()
                    .delimiter(b'\t')
                    .has_headers(true)
                    .from_path(coverage_table_path)?;
                let headers = reader.headers()?;
                // skip first three columns, then take every other column starting from fourth
                for header in headers.iter().skip(3).step_by(2) {
                    if header.contains("/") {
                        // when performing read mapping, coverm sets the column name to
                        // {reference}/{sample}.bam
                        let mut sample_name = header.split("/").last().unwrap().to_string();
                        sample_name = sample_name.replace(".bam", "");
                        previous_sample_names.insert(sample_name);
                    } else {
                        // when using BAM files, coverm sets the column name to
                        // {sample}
                        previous_sample_names.insert(header.to_string());
                    }
                }
                Ok(Some(previous_sample_names))
            },
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
    pub fn new(m: &clap::ArgMatches) -> Result<Self> {
        let mut read1: Option<Vec<_>> = None;
        let mut read2: Option<Vec<_>> = None;
        let mut interleaved: Option<Vec<_>> = None;
        let mut unpaired: Option<Vec<_>> = None;
        let mut long_reads: Option<Vec<_>> = None;
        let mut short_read_bams: Option<Vec<_>> = None;
        let mut long_read_bams: Option<Vec<_>> = None;

        if m.contains_id("read1") {
            let inner_read1: Vec<_> = m
                .get_many::<String>("read1")
                .unwrap()
                .map(|s| s.clone())
                .collect();
            let inner_read2: Vec<_> = m
                .get_many::<String>("read2")
                .unwrap()
                .map(|s| s.clone())
                .collect();
            if inner_read1.len() != inner_read2.len() {
                return Err(anyhow!(
                    "When specifying paired reads with the -1 and -2 flags, \
                        there must be equal numbers specified. Instead found \
                        {} and {} respectively",
                    inner_read1.len(),
                    inner_read2.len()
                ));
            }
            read1 = Some(inner_read1);
            read2 = Some(inner_read2);
        }

        // Parse --coupled
        if m.contains_id("coupled") {
            let coupled: Vec<_> = m
                .get_many::<String>("coupled")
                .unwrap()
                .map(|s| s.clone())
                .collect();
            if coupled.len() % 2 != 0 {
                return Err(anyhow!(
                    "The --coupled flag must be set with pairs of read \
                     sets, but an odd number ({}) was specified",
                    coupled.len()
                ));
            }
            let mut i = 0;
            let mut inner_read1 = Vec::with_capacity(coupled.len() / 2);
            let mut inner_read2 = Vec::with_capacity(coupled.len() / 2);
            while i < coupled.len() {
                inner_read1.push(coupled[i].clone());
                inner_read2.push(coupled[i + 1].clone());
                i += 2;
            }

            if read1.is_some() {
                read1.as_mut().unwrap().extend(inner_read1);
                read2.as_mut().unwrap().extend(inner_read2);
            } else {
                read1 = Some(inner_read1);
                read2 = Some(inner_read2);
            }
        }

        if m.contains_id("interleaved") {
            interleaved = Some(m
                .get_many::<String>("interleaved")
                .unwrap()
                .map(|s| s.clone())
                .collect());
        }
        if m.contains_id("single") {
            unpaired = Some(m
                .get_many::<String>("single")
                .unwrap()
                .map(|s| s.clone())
                .collect());
        }

        if m.contains_id("longreads") {
            long_reads = Some(m
                .get_many::<String>("longreads")
                .unwrap()
                .map(|s| s.clone())
                .collect());
        }

        if m.contains_id("longread-bam-files") {
            long_read_bams = Some(m
                .get_many::<String>("longread-bam-files")
                .unwrap()
                .map(|s| s.clone())
                .collect());
        }

        if m.contains_id("bam-files") {
            short_read_bams = Some(m
                .get_many::<String>("bam-files")
                .unwrap()
                .map(|s| s.clone())
                .collect());
        }

        Ok(
            Self {
                forward_read_paths: read1,
                reverse_read_paths: read2,
                interleaved_read_paths: interleaved,
                unpaired_read_paths: unpaired,
                long_read_paths: long_reads,
                short_read_bam_paths: short_read_bams,
                long_read_bam_paths: long_read_bams,
            }
        )
    }

    /// returns a ReadCollection that only contains the
    /// sample names specified in the sample_names_to_map
    /// subset to only non-longread samples and no BAM files
    pub fn subset_short_reads(&self, sample_names_to_map: &HashSet<&str>) -> Self {
        let mut read1: Option<Vec<_>> = None;
        let mut read2: Option<Vec<_>> = None;
        let mut interleaved: Option<Vec<_>> = None;
        let mut unpaired: Option<Vec<_>> = None;

        if let Some(read1_paths) = &self.forward_read_paths {
            let mut inner_read1 = Vec::new();
            let mut inner_read2 = Vec::new();
            for (read1_path, read2_path) in read1_paths.iter().zip(self.reverse_read_paths.as_ref().unwrap()) {
                let read1_path = Path::new(read1_path);
                let read2_path = Path::new(read2_path);
                let sample_name = read1_path.file_name().unwrap().to_str().unwrap();
                if sample_names_to_map.contains(sample_name) {
                    inner_read1.push(read1_path.to_str().unwrap().to_string());
                    inner_read2.push(read2_path.to_str().unwrap().to_string());
                }
            }
            read1 = Some(inner_read1);
            read2 = Some(inner_read2);
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
            interleaved = Some(inner_interleaved);
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
            unpaired = Some(inner_unpaired);
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
            let mut inner_long_reads = Vec::new();
            for long_read_path in long_read_paths {
                let long_read_path = Path::new(long_read_path);
                let sample_name = long_read_path.file_name().unwrap().to_str().unwrap();
                if sample_names_to_map.contains(sample_name) {
                    inner_long_reads.push(long_read_path.to_str().unwrap().to_string());
                }
            }
            long_reads = Some(inner_long_reads);
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
            short_read_bams = Some(inner_short_read_bams);
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
            long_read_bams = Some(inner_long_read_bams);
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
            coverm_command.arg("--unpaired");
            for path in unpaired {
                coverm_command.arg(path);
            }
            return;
        }

        if let Some(long_reads) = &self.long_read_paths {
            coverm_command.arg("--unpaired");
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
            return;
        }
    }
}

const EPSILON: f64 = 0.0000001;
const MIN_VAR_EPSILON: f64 = 1e-4;
const MIN_VAR: f64 = 1.0;
pub struct MetabatDistance;

impl MetabatDistance {
    pub fn distance<D: Dimension>(coverage_array1: ArrayView<f64, D>, coverage_array2: ArrayView<f64, D>) -> f64 {
        
        let n_samples = coverage_array1.len() / 2;
        let mut mb_vec = Vec::with_capacity(n_samples);

        let a_means = coverage_array1.iter().step_by(2);
        let b_means = coverage_array2.iter().step_by(2);
        let a_vars = coverage_array1.iter().skip(1).step_by(2);
        let b_vars = coverage_array2.iter().skip(1).step_by(2);

        let mut both_present = Vec::with_capacity(n_samples);
        for (sample_index, (a_mean, b_mean, a_var, b_var)) in izip!(a_means, b_means, a_vars, b_vars).enumerate() {
            let a_mean = a_mean + &EPSILON;
            let b_mean = b_mean + &EPSILON;
            let mut a_var = a_var + &EPSILON;
            let mut b_var = b_var + &EPSILON;
            if a_var < MIN_VAR {
                a_var = MIN_VAR;
            }
            if b_var < MIN_VAR {
                b_var = MIN_VAR;
            }


            if a_mean > EPSILON && b_mean > EPSILON {
                both_present.push(sample_index);
            }

            let mut k1;
            let mut k2;

            if (a_mean > EPSILON || b_mean > EPSILON) && a_mean != b_mean {
                if (a_var - b_var).abs() < MIN_VAR_EPSILON {
                    let tmp = ((a_mean + b_mean) / 2.0) as f64;
                    k1 = tmp;
                    k2 = tmp;
                } else {
                    let tmp = (a_var * b_var * ((a_mean - b_mean) * (a_mean - b_mean) - 2.0 * (a_var - b_var) * ((b_var / a_var).sqrt()).ln())).sqrt();
                    k1 = ((tmp - a_mean * b_var + b_mean * a_var) / (a_var - b_var)) as f64;
                    k2 = ((tmp + a_mean * b_var - b_mean * a_var) / (b_var - a_var)) as f64;
                }

                if k1 > k2 {
                    std::mem::swap(&mut k1, &mut k2);
                }

                let p1;
                let p2;

                // may not need to be calculate square root of variances here?
                if a_var > b_var {
                    p1 = Normal::new(b_mean as f64, b_var.sqrt() as f64).unwrap();
                    p2 = Normal::new(a_mean as f64, a_var.sqrt() as f64).unwrap();
                } else {
                    p1 = Normal::new(a_mean as f64, a_var.sqrt() as f64).unwrap();
                    p2 = Normal::new(b_mean as f64, b_var.sqrt() as f64).unwrap();
                }

                if (k1 - k2).abs() < EPSILON as f64 {
                    let d = p1.cdf(k1) - p2.cdf(k1);
                    mb_vec.push(d as f64);
                } else {
                    let d = (p1.cdf(k2) - p1.cdf(k1) + p2.cdf(k1) - p2.cdf(k2)).abs();
                    mb_vec.push(d as f64);
                }
            }
        }

        if mb_vec.len() > 0 {
            // ln mean of mb_vec
            let length = mb_vec.len();
            let mut d = (mb_vec.into_iter().map(|x| x.ln()).sum::<f64>() / length as f64).exp();
            if d.is_nan() {
                d = 1.0;
            }
            return d
        }
        return 1.0
    }
}

// impl Distance<D: Clone + Send + Sync> for MetabatDistance {
//     fn 
// }

impl Distance<f64> for MetabatDistance {
    fn eval(&self, coverage_array1: &[f64], coverage_array2: &[f64]) -> f32 {
        let n_samples = coverage_array1.len() / 2;
        let mut mb_vec = Vec::with_capacity(n_samples);

        let a_means = coverage_array1.iter().step_by(2);
        let b_means = coverage_array2.iter().step_by(2);
        let a_vars = coverage_array1.iter().skip(1).step_by(2);
        let b_vars = coverage_array2.iter().skip(1).step_by(2);

        let mut both_present = Vec::with_capacity(n_samples);
        for (sample_index, (a_mean, b_mean, a_var, b_var)) in izip!(a_means, b_means, a_vars, b_vars).enumerate() {
            let a_mean = a_mean + &EPSILON;
            let b_mean = b_mean + &EPSILON;
            let mut a_var = a_var + &EPSILON;
            let mut b_var = b_var + &EPSILON;
            if a_var < MIN_VAR {
                a_var = MIN_VAR;
            }
            if b_var < MIN_VAR {
                b_var = MIN_VAR;
            }


            if a_mean > EPSILON && b_mean > EPSILON {
                both_present.push(sample_index);
            }

            let mut k1;
            let mut k2;

            if (a_mean > EPSILON || b_mean > EPSILON) && a_mean != b_mean {
                if (a_var - b_var).abs() < MIN_VAR_EPSILON {
                    let tmp = ((a_mean + b_mean) / 2.0) as f64;
                    k1 = tmp;
                    k2 = tmp;
                } else {
                    let tmp = (a_var * b_var * ((a_mean - b_mean) * (a_mean - b_mean) - 2.0 * (a_var - b_var) * ((b_var / a_var).sqrt()).ln())).sqrt();
                    k1 = ((tmp - a_mean * b_var + b_mean * a_var) / (a_var - b_var)) as f64;
                    k2 = ((tmp + a_mean * b_var - b_mean * a_var) / (b_var - a_var)) as f64;
                }

                if k1 > k2 {
                    std::mem::swap(&mut k1, &mut k2);
                }

                let p1;
                let p2;

                // may not need to be calculate square root of variances here?
                if a_var > b_var {
                    p1 = Normal::new(b_mean as f64, b_var.sqrt() as f64).unwrap();
                    p2 = Normal::new(a_mean as f64, a_var.sqrt() as f64).unwrap();
                } else {
                    p1 = Normal::new(a_mean as f64, a_var.sqrt() as f64).unwrap();
                    p2 = Normal::new(b_mean as f64, b_var.sqrt() as f64).unwrap();
                }

                if (k1 - k2).abs() < EPSILON as f64 {
                    let d = p1.cdf(k1) - p2.cdf(k1);
                    mb_vec.push(d as f64);
                } else {
                    let d = (p1.cdf(k2) - p1.cdf(k1) + p2.cdf(k1) - p2.cdf(k2)).abs();
                    mb_vec.push(d as f64);
                }
            }
        }

        if mb_vec.len() > 0 {
            // ln mean of mb_vec
            let length = mb_vec.len();
            let d = (mb_vec.into_iter().map(|x| x.ln()).sum::<f64>() / length as f64).exp();
            return d as f32
        }
        return 1.0
    }
}