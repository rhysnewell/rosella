pub mod booster;
pub mod tables;

use std::collections::HashMap;
use std::io::{BufRead, BufWriter, Read, Write};
use std::path::Path;

use anyhow::Result;
use log::info;

use crate::external::diamond_engine::DiamondEngine;
use crate::markers::orfs;
use booster::Booster;
use tables::{METADATA, Tables};

const COMPLETENESS_GZ: &[u8] = include_bytes!("../../data/checkm2_completeness.gbm.gz");
const CONTAMINATION_GZ: &[u8] = include_bytes!("../../data/checkm2_contamination.gbm.gz");
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

/// Gene families per contig, summed into a bin on demand. The trained columns are counts and
/// residue totals, so a candidate's vector is the sum of its contigs and the assembly is
/// annotated once however many candidates the search proposes.
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

impl ContigQuality {
    pub fn annotate(
        assembly: &str,
        names: &[String],
        threads: usize,
        database: &Path,
    ) -> Result<Self> {
        let tables = Tables::load()?;
        let engine = DiamondEngine::new(database, threads)?;
        let index = names
            .iter()
            .enumerate()
            .map(|(position, name)| (name.as_str(), position))
            .collect::<HashMap<_, _>>();

        info!("Calling genes over the assembly.");
        let contigs = orfs::read_wanted(assembly, &index)?;
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
        engine.best_hits(&proteins, &table)?;

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

        Ok(Self {
            tables,
            completeness: Booster::parse(&inflate(COMPLETENESS_GZ)?)?,
            contamination: Booster::parse(&inflate(CONTAMINATION_GZ)?)?,
            metadata,
            hits,
        })
    }

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
