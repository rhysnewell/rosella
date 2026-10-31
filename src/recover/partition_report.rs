use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::Result;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::{Graph, row_of};

pub const SAMPLE_PASS: &str = "sample";
pub const PULL: &str = "pull";
const KEPT_PASS: &str = "derived";

pub struct PartitionReport {
    contigs: Vec<usize>,
    pass: String,
    headers: Vec<String>,
    columns: Vec<Vec<Option<usize>>>,
}

impl PartitionReport {
    pub fn new(contigs: &[usize]) -> Self {
        Self {
            contigs: contigs.to_vec(),
            pass: SAMPLE_PASS.to_string(),
            headers: Vec::new(),
            columns: Vec::new(),
        }
    }

    pub fn pass(&mut self, pass: &str) {
        self.pass = pass.to_string();
    }

    // Readers take the derived block as the partition the run used, so only the kept pass
    // carries that name and every refused pass keeps its own.
    pub fn keep(&mut self, pass: &str) {
        let prefix = format!("{pass}_");
        for header in &mut self.headers {
            if let Some(rest) = header.strip_prefix(&prefix) {
                *header = format!("{KEPT_PASS}_{rest}");
            }
        }
    }

    pub fn ladder(&mut self, ladder: &[Partitioning]) {
        let mut rungs: HashMap<(&str, u64), usize> = HashMap::new();
        for held in ladder {
            let rung = rungs.entry((held.arm.name(), held.seed)).or_default();
            let name = format!("{}_s{}_r{rung}", held.arm.name(), held.seed);
            *rung += 1;
            self.add(&name, held);
        }
    }

    pub fn chosen(&mut self, arms: &[Partitioning]) {
        let mut kept: HashMap<(&str, u64), usize> = HashMap::new();
        for held in arms {
            let at = kept.entry((held.arm.name(), held.seed)).or_default();
            let name = format!("{}_s{}_{at}_chosen", held.arm.name(), held.seed);
            *at += 1;
            self.add(&name, held);
        }
    }

    pub fn add(&mut self, name: &str, held: &Partitioning) {
        self.headers.push(format!("{}_{name}", self.pass));
        self.columns.push(labels(held, self.contigs.len()));
    }

    pub fn write(&self, path: &Path, names: &[String], lengths: &[usize]) -> Result<()> {
        let mut out = BufWriter::new(File::create(path)?);
        writeln!(out, "contig\tlength\t{}", self.headers.join("\t"))?;
        for (at, contig) in self.contigs.iter().enumerate() {
            write!(out, "{}\t{}", names[*contig], lengths[*contig])?;
            for column in &self.columns {
                match column[at] {
                    Some(label) => write!(out, "\t{label}")?,
                    None => write!(out, "\t-")?,
                }
            }
            writeln!(out)?;
        }
        out.flush()?;
        Ok(())
    }
}

fn labels(held: &Partitioning, positions: usize) -> Vec<Option<usize>> {
    let mut labels = vec![None; positions];
    for (label, members) in &held.cluster_map {
        for at in members {
            labels[*at] = Some(*label);
        }
    }
    labels
}

pub fn pass_file(path: &Path, pass: &str) -> PathBuf {
    let mut name = path.file_name().unwrap_or_default().to_os_string();
    name.push(format!(".{pass}"));
    path.with_file_name(name)
}

// Each pass writes its own report, so the one left at the path is the kept pass.
pub fn keep_pass(path: &Path, kept: &str, passes: &[String]) -> Result<()> {
    for pass in passes {
        let file = pass_file(path, pass);
        match pass == kept {
            true => std::fs::rename(&file, path)?,
            false => std::fs::remove_file(&file)?,
        }
    }
    Ok(())
}

// Leiden reads the bp-weighted graph, so both sums are written to tell a pull the lengths make
// from one the distances make.
pub fn write_pull(
    path: &Path,
    graph: &Graph,
    held: &Partitioning,
    contigs: &[usize],
    names: &[String],
    lengths: &[usize],
) -> Result<()> {
    let community = labels(held, contigs.len());
    let label = |held: Option<usize>| held.map_or("-".to_string(), |label| label.to_string());
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "contig\tcommunity\ttoward\tweight\tbp_weight")?;
    let mut pull: BTreeMap<Option<usize>, (f64, f64)> = BTreeMap::new();
    for (at, contig) in contigs.iter().enumerate() {
        pull.clear();
        let (targets, weights) = row_of(graph, at);
        for (target, weight) in targets.iter().zip(weights) {
            let target = *target as usize;
            if target == at {
                continue;
            }
            let bases = (lengths[*contig] as f64 * lengths[contigs[target]] as f64).sqrt();
            let sums = pull.entry(community[target]).or_default();
            sums.0 += *weight as f64;
            sums.1 += *weight as f64 * bases;
        }
        for (toward, (weight, bp_weight)) in &pull {
            writeln!(
                out,
                "{}\t{}\t{}\t{weight:.6}\t{bp_weight:.1}",
                names[*contig],
                label(community[at]),
                label(*toward)
            )?;
        }
    }
    out.flush()?;
    Ok(())
}
