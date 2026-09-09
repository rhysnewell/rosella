pub mod booster;
pub mod cache;
pub mod orfs;
pub mod tables;

use std::collections::HashMap;
use std::io::{BufRead, BufWriter, Read, Write};
use std::path::Path;

use anyhow::Result;
use log::{info, warn};

use crate::external::diamond_engine::{DiamondEngine, Sensitivity};
use booster::Booster;
use tables::{METADATA, Tables};

const COMPLETENESS_GZ: &[u8] = include_bytes!("../../data/completeness.gbm.gz");
const CONTAMINATION_GZ: &[u8] = include_bytes!("../../data/contamination.gbm.gz");
const RESIDUES: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

#[derive(Debug, Clone, Copy, Default)]
pub struct Quality {
    pub completeness: f64,
    pub contamination: f64,
}

impl Quality {
    /// The standard genome score, so one number orders candidates that trade the two.
    pub fn score(&self) -> f64 {
        self.completeness - 5.0 * self.contamination
    }
}

pub trait Scorer: Sync {
    fn score(&self, contigs: &[usize]) -> Quality;
}

impl Scorer for ContigQuality {
    fn score(&self, contigs: &[usize]) -> Quality {
        ContigQuality::score(self, contigs)
    }
}

/// Every trained column is a per contig sum, so the assembly is annotated once however many
/// candidates the search proposes.
pub struct ContigQuality {
    tables: Tables,
    completeness: Booster,
    contamination: Booster,
    metadata: Vec<[u32; METADATA]>,
    hits: Vec<Vec<(u32, u32)>>,
}

fn inflate(compressed: &[u8]) -> Result<String> {
    let mut text = String::new();
    flate2::read::GzDecoder::new(compressed).read_to_string(&mut text)?;
    Ok(text)
}

fn residue_column(residue: u8) -> Option<usize> {
    RESIDUES.iter().position(|wanted| *wanted == residue)
}

fn write_proteins(orfs: &[orfs::Orf], target: &Path) -> Result<()> {
    let mut sink = BufWriter::new(std::fs::File::create(target)?);
    for (position, orf) in orfs.iter().enumerate() {
        if orf.protein.is_empty() {
            continue;
        }
        writeln!(sink, ">{position}\n{}", orf.protein)?;
    }
    sink.flush()?;
    Ok(())
}

fn restore(path: &std::path::Path) -> Option<(Vec<String>, cache::Annotation)> {
    if !path.is_file() {
        return None;
    }
    match cache::read(path) {
        Ok(held) => {
            info!("Read gene families from {}.", path.display());
            Some(held)
        }
        Err(error) => {
            warn!("Ignoring {}: {error}", path.display());
            None
        }
    }
}

fn search(
    assembly: &str,
    min_contig_size: usize,
    threads: usize,
    sensitivity: Sensitivity,
    database: &Path,
    tables: &Tables,
) -> Result<(Vec<String>, cache::Annotation)> {
    let engine = DiamondEngine::new(database, threads, sensitivity)?;

    info!("Calling genes over the assembly.");
    let (names, contigs) = orfs::read_over(assembly, min_contig_size)?;
    let orfs = orfs::call(&contigs, threads)?;

    let mut metadata = vec![[0u32; METADATA]; names.len()];
    for orf in &orfs {
        let row = &mut metadata[orf.contig];
        for residue in orf.protein.bytes() {
            if let Some(column) = residue_column(residue) {
                row[column] += 1;
            }
        }
        // The reference pipeline reads a protein file that still carries the stop.
        row[20] += orf.protein.len() as u32 + u32::from(!orf.partial);
        row[21] += 1;
    }

    let workspace = tempfile::tempdir()?;
    let proteins = workspace.path().join("proteins.faa");
    let table = workspace.path().join("hits.tsv");
    write_proteins(&orfs, &proteins)?;
    info!("Searching {} proteins for gene families.", orfs.len());
    engine.best_hits(&proteins, &table, workspace.path())?;

    let mut counts = vec![HashMap::<u32, u32>::new(); names.len()];
    let reader = std::io::BufReader::new(std::fs::File::open(&table)?);
    for line in reader.lines() {
        let line = line?;
        let mut fields = line.split('\t');
        let (Some(query), Some(subject)) = (fields.next(), fields.next()) else {
            continue;
        };
        let Ok(position) = query.parse::<usize>() else {
            continue;
        };
        let Some(family) = subject.split('~').nth(1) else {
            continue;
        };
        if let Some(column) = tables.kos.get(family) {
            *counts[orfs[position].contig].entry(*column).or_default() += 1;
        }
    }

    let hits = counts
        .into_iter()
        .map(|held| {
            let mut held = held.into_iter().collect::<Vec<_>>();
            held.sort_unstable();
            held
        })
        .collect();

    Ok((names, cache::Annotation { metadata, hits }))
}

/// The search runs before the contigs a run keeps are known, so it annotates every contig over
/// the length floor and the run takes the rows it needs out afterwards.
pub struct Annotated {
    tables: Tables,
    names: Vec<String>,
    annotation: cache::Annotation,
}

impl Annotated {
    pub fn build(
        assembly: &str,
        min_contig_size: usize,
        threads: usize,
        sensitivity: Sensitivity,
        database: &Path,
        cache_directory: Option<&Path>,
    ) -> Result<Self> {
        let tables = Tables::load()?;
        let stored = cache_directory
            .map(|home| cache::path_for(home, assembly, database, sensitivity.tag()));
        let (names, annotation) = match stored.as_deref().and_then(restore) {
            Some(held) => held,
            None => {
                let (names, annotation) = search(
                    assembly,
                    min_contig_size,
                    threads,
                    sensitivity,
                    database,
                    &tables,
                )?;
                if let Some(path) = stored.as_deref()
                    && let Err(error) = cache::write(path, &names, &annotation)
                {
                    warn!("Could not cache the gene families: {error}");
                }
                (names, annotation)
            }
        };
        Ok(Self {
            tables,
            names,
            annotation,
        })
    }

    pub fn select(self, names: &[String]) -> Result<ContigQuality> {
        let annotation = cache::select(names, &self.names, self.annotation)?;
        Ok(ContigQuality {
            tables: self.tables,
            completeness: Booster::parse(&inflate(COMPLETENESS_GZ)?)?,
            contamination: Booster::parse(&inflate(CONTAMINATION_GZ)?)?,
            metadata: annotation.metadata,
            hits: annotation.hits,
        })
    }
}

impl ContigQuality {
    pub fn write_report(
        &self,
        bins: &std::collections::BTreeMap<usize, Vec<usize>>,
        lengths: &[usize],
        path: &Path,
    ) -> Result<()> {
        let mut sink = BufWriter::new(std::fs::File::create(path)?);
        writeln!(sink, "bin\tcontigs\tbp\tcompleteness\tcontamination")?;
        for (bin, contigs) in bins {
            let held = self.score(contigs);
            let bp = contigs.iter().map(|contig| lengths[*contig]).sum::<usize>();
            writeln!(
                sink,
                "rosella_bin_{bin}\t{}\t{bp}\t{:.2}\t{:.2}",
                contigs.len(),
                held.completeness,
                held.contamination
            )?;
        }
        sink.flush()?;
        Ok(())
    }

    /// A partner that brings no family the bin lacks cannot raise its completeness, which is
    /// most pairs, so this keeps the search off the boosters.
    pub fn families(&self, contigs: &[usize]) -> std::collections::HashSet<u32> {
        contigs
            .iter()
            .flat_map(|contig| self.hits[*contig].iter().map(|(column, _)| *column))
            .collect()
    }

    pub fn score(&self, contigs: &[usize]) -> Quality {
        let mut metadata = [0.0f64; METADATA];
        let mut genes = vec![0.0f64; self.tables.gene_count];
        for contig in contigs {
            for (column, count) in self.metadata[*contig].iter().enumerate() {
                metadata[column] += *count as f64;
            }
            for (column, count) in &self.hits[*contig] {
                genes[*column as usize] += *count as f64;
            }
        }

        let mut vector = Vec::with_capacity(self.tables.width());
        self.tables.fill(&metadata, &genes, &mut vector);
        Quality {
            completeness: self.completeness.predict(&vector).clamp(0.0, 100.0),
            contamination: self.contamination.predict(&vector).max(0.0),
        }
    }
}
