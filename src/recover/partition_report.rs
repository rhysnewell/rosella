use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::clustering::clusterer::Partitioning;

pub struct PartitionReport {
    contigs: Vec<usize>,
    pass: &'static str,
    headers: Vec<String>,
    columns: Vec<Vec<Option<usize>>>,
}

impl PartitionReport {
    pub fn new(contigs: &[usize]) -> Self {
        Self {
            contigs: contigs.to_vec(),
            pass: "sample",
            headers: Vec::new(),
            columns: Vec::new(),
        }
    }

    pub fn pass(&mut self, pass: &'static str) {
        self.pass = pass;
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
        for held in arms {
            self.add(&format!("{}_s{}_chosen", held.arm.name(), held.seed), held);
        }
    }

    pub fn add(&mut self, name: &str, held: &Partitioning) {
        let mut labels = vec![None; self.contigs.len()];
        for (label, members) in &held.cluster_map {
            for at in members {
                labels[*at] = Some(*label);
            }
        }
        self.headers.push(format!("{}_{name}", self.pass));
        self.columns.push(labels);
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
