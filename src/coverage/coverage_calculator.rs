use std::collections::HashSet;
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use log::info;

use super::bam_depth::{self, Sample};
use super::coverage_table::{CoverageTable, bam_stem};
use super::depth_file::{self, CONTIGS};
use crate::cli::{CoverageSource, ReadArgs};
use crate::external::coverm_engine::Mapper;

pub struct CoverageInputs<'a> {
    pub assembly: &'a str,
    pub output_directory: &'a str,
    pub depth: &'a Path,
    pub threads: usize,
    pub reads: &'a ReadArgs,
}

pub struct Coverage {
    pub table: CoverageTable,
    pub file: String,
    pub windows: Option<PathBuf>,
}

pub fn samples(inputs: &CoverageInputs) -> Result<Vec<Sample>> {
    let collection = ReadCollection::new(&inputs.reads.coverage)?;
    let mut samples = Mapper::new(inputs).bams(&collection.reads)?;
    samples.extend(collection.bams.iter().map(|path| Sample {
        name: bam_stem(path).to_string(),
        path: path.into(),
    }));
    let mut names = HashSet::new();
    if let Some(sample) = samples.iter().find(|sample| !names.insert(&sample.name)) {
        bail!(
            "two inputs name the sample {}, rename one file",
            sample.name
        );
    }
    Ok(samples)
}

// The file returned holds exactly that table, so rows left out can be read back later.
pub fn calculate_coverage(inputs: &CoverageInputs) -> Result<Coverage> {
    std::fs::create_dir_all(inputs.output_directory)?;
    let source = &inputs.reads.coverage;
    let samples = samples(inputs)?;
    let read = match (samples.is_empty(), &source.depth_directory) {
        (false, _) => Some(bam_depth::read(&samples, inputs.depth, inputs.threads)?),
        (true, Some(_)) if source.coverage_file.is_none() => Some(depth_file::kept(inputs.depth)?),
        (true, _) => None,
    };
    let (table, file) = match (read, &source.coverage_file) {
        (None, Some(held)) => (CoverageTable::from_any_file(held)?, held.clone()),
        (None, None) => bail!("no coverage table, reads or BAMs were given"),
        (Some(mut table), held) => {
            if let Some(held) = held {
                let mut merged = CoverageTable::from_any_file(held)?;
                merged.merge(table)?;
                table = merged;
            }
            let file = format!("{}/coverage.tsv", inputs.output_directory);
            table.write(&file)?;
            (table, file)
        }
    };
    let windows = windows_for(inputs.depth, &table.sample_names);
    Ok(Coverage {
        table,
        file,
        windows,
    })
}

fn windows_for(depth: &Path, samples: &[String]) -> Option<PathBuf> {
    let missing = samples
        .iter()
        .filter(|name| !depth_file::path(depth, name).is_file())
        .count();
    if missing == 0 && depth.join(CONTIGS).is_file() {
        return Some(depth.to_path_buf());
    }
    if missing < samples.len() {
        info!(
            "Only {} of {} samples have window depth in {}, so it goes unused.",
            samples.len() - missing,
            samples.len(),
            depth.display()
        );
    }
    None
}

pub enum Reads {
    Paired(String, String),
    Interleaved(String),
    Single(String),
    Long(String),
}

impl Reads {
    pub fn first(&self) -> &str {
        match self {
            Self::Paired(first, _)
            | Self::Interleaved(first)
            | Self::Single(first)
            | Self::Long(first) => first,
        }
    }

    pub fn name(&self) -> &str {
        let name = bam_stem(self.first());
        let name = [".gz", ".bz2", ".xz", ".zst"]
            .iter()
            .find_map(|suffix| name.strip_suffix(suffix))
            .unwrap_or(name);
        [".fastq", ".fq", ".fasta", ".fa", ".fna"]
            .iter()
            .find_map(|suffix| name.strip_suffix(suffix))
            .unwrap_or(name)
    }

    pub fn is_long(&self) -> bool {
        matches!(self, Self::Long(_))
    }
}

pub struct ReadCollection {
    pub reads: Vec<Reads>,
    pub bams: Vec<String>,
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
        let paired = source
            .read1
            .iter()
            .zip(&source.read2)
            .chain(source.coupled.chunks(2).map(|pair| (&pair[0], &pair[1])))
            .map(|(forward, reverse)| Reads::Paired(forward.clone(), reverse.clone()));
        let reads = paired
            .chain(source.interleaved.iter().cloned().map(Reads::Interleaved))
            .chain(source.single.iter().cloned().map(Reads::Single))
            .chain(source.longreads.iter().cloned().map(Reads::Long))
            .collect();
        let bams = source
            .bam_files
            .iter()
            .chain(&source.longread_bam_files)
            .cloned()
            .collect();
        Ok(Self { reads, bams })
    }
}
