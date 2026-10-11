use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use ndarray::Array2;

use super::coverage_table::CoverageTable;
use crate::learned::windows::WINDOW;

pub const CONTIGS: &str = "contigs.tsv";
pub const TABLE: &str = "coverage.tsv.gz";
const SUFFIX: &str = ".depth.zst";
const PER_CONTIG: usize = 8 + 8 + 4;
const LEVEL: i32 = 3;

#[derive(PartialEq)]
pub struct Contigs {
    pub names: Vec<String>,
    pub lengths: Vec<usize>,
}

pub struct Depth {
    pub means: Vec<f64>,
    pub variances: Vec<f64>,
    pub reads: Vec<u32>,
}

pub fn path(directory: &Path, name: &str) -> PathBuf {
    directory.join(format!("{name}{SUFFIX}"))
}

impl Contigs {
    pub fn read(directory: &Path) -> Result<Self> {
        let path = directory.join(CONTIGS);
        let (mut names, mut lengths) = (Vec::new(), Vec::new());
        let file = File::open(&path).with_context(|| format!("opening {}", path.display()))?;
        for line in BufReader::new(file).lines().skip(1) {
            let line = line?;
            let (name, length) = line.split_once('\t').context("a contig has no length")?;
            names.push(name.to_string());
            lengths.push(length.parse()?);
        }
        Ok(Self { names, lengths })
    }

    pub fn keep(&self, directory: &Path) -> Result<()> {
        if directory.join(CONTIGS).is_file() {
            if *self != Self::read(directory)? {
                bail!(
                    "{} holds the depth of other contigs than these BAMs were mapped to, so \
                     remove it or name another --depth-directory",
                    directory.display()
                );
            }
            return Ok(());
        }
        let rows = self
            .names
            .iter()
            .zip(&self.lengths)
            .map(|(name, length)| format!("{name}\t{length}\n"))
            .collect::<String>();
        crate::report_sink::write_atomically(&directory.join(CONTIGS), |mut file| {
            Ok(file.write_all(format!("contig\tlength\n{rows}").as_bytes())?)
        })
    }

    pub fn windows(&self) -> usize {
        self.lengths
            .iter()
            .map(|length| length.div_ceil(WINDOW))
            .sum()
    }
}

// The per-contig columns come before the windows, so a table is read without decoding them.
pub fn write(path: &Path, depth: &Depth, windows: &[u8]) -> Result<()> {
    crate::report_sink::write_atomically(path, |file| {
        let mut encoder = zstd::Encoder::new(file, LEVEL)?;
        let columns = depth.means.iter().chain(&depth.variances);
        for value in columns {
            encoder.write_all(&value.to_le_bytes())?;
        }
        for reads in &depth.reads {
            encoder.write_all(&reads.to_le_bytes())?;
        }
        encoder.write_all(windows)?;
        encoder.finish()?;
        Ok(())
    })
}

pub fn read(path: &Path, contigs: usize) -> Result<Depth> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut bytes = vec![0; PER_CONTIG * contigs];
    zstd::Decoder::new(file)?
        .read_exact(&mut bytes)
        .with_context(|| format!("{} holds fewer than {contigs} contigs", path.display()))?;
    let (floats, reads) = bytes.split_at(16 * contigs);
    let floats = floats
        .chunks_exact(8)
        .map(|value| f64::from_le_bytes(value.try_into().expect("eight bytes")))
        .collect::<Vec<_>>();
    let (means, variances) = floats.split_at(contigs);
    Ok(Depth {
        means: means.to_vec(),
        variances: variances.to_vec(),
        reads: reads
            .chunks_exact(4)
            .map(|value| u32::from_le_bytes(value.try_into().expect("four bytes")))
            .collect(),
    })
}

pub fn windows(path: &Path, contigs: &Contigs) -> Result<Vec<u8>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut bytes = Vec::new();
    zstd::Decoder::new(file)?.read_to_end(&mut bytes)?;
    let (head, windows) = (PER_CONTIG * contigs.names.len(), 2 * contigs.windows());
    if bytes.len() != head + windows {
        bail!(
            "{} holds {} bytes where {} contigs and {} windows need {}",
            path.display(),
            bytes.len(),
            contigs.names.len(),
            windows / 2,
            head + windows
        );
    }
    bytes.drain(..head);
    Ok(bytes)
}

pub fn table(directory: &Path, contigs: &Contigs, name: &str) -> Result<CoverageTable> {
    let depth = read(&path(directory, name), contigs.names.len())?;
    let mut table = Array2::zeros((contigs.names.len(), 2));
    for (row, (mean, variance)) in depth.means.iter().zip(&depth.variances).enumerate() {
        table[[row, 0]] = *mean;
        table[[row, 1]] = *variance;
    }
    Ok(CoverageTable::new(
        contigs.names.clone(),
        contigs.lengths.clone(),
        vec![name.to_string()],
        table,
    ))
}

pub fn merged(directory: &Path, contigs: &Contigs, names: &[String]) -> Result<CoverageTable> {
    let mut tables = names.iter().map(|name| table(directory, contigs, name));
    let mut merged = tables
        .next()
        .with_context(|| format!("{} holds no sample depth", directory.display()))??;
    for next in tables {
        merged.merge(next?)?;
    }
    Ok(merged)
}

// Samples are taken in name order, so a directory filled one sample at a time reads the same
// however it was filled.
pub fn kept(directory: &Path) -> Result<CoverageTable> {
    let mut names = std::fs::read_dir(directory)
        .with_context(|| format!("reading {}", directory.display()))?
        .map(|entry| Ok(entry?.file_name()))
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .filter_map(|name| name.into_string().ok())
        // macOS leaves a hidden ._ twin beside each file on drives without extended attributes.
        .filter(|name| !name.starts_with('.'))
        .filter_map(|name| name.strip_suffix(SUFFIX).map(str::to_string))
        .collect::<Vec<_>>();
    names.sort();
    merged(directory, &Contigs::read(directory)?, &names)
}
