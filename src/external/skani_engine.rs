use std::{
    collections::HashMap,
    io::{BufWriter, Write},
    path::Path,
    process::Command,
};

use anyhow::Result;
use log::{debug, info};
use needletail::parse_fastx_file;

use crate::homology::Pair;

pub struct SkaniEngine<'a> {
    assembly: &'a str,
    threads: usize,
    min_contig_length: usize,
}

impl<'a> SkaniEngine<'a> {
    pub fn new(assembly: &'a str, threads: usize, min_contig_length: usize) -> Self {
        Self {
            assembly,
            threads,
            min_contig_length,
        }
    }

    pub fn is_installed() -> bool {
        Command::new("skani")
            .arg("--version")
            .output()
            .is_ok_and(|output| output.status.success())
    }

    pub fn pairs(&self, names: &[String], lengths: &[usize]) -> Result<Vec<Pair>> {
        let wanted = names
            .iter()
            .zip(lengths)
            .filter(|(_, length)| **length >= self.min_contig_length)
            .map(|(name, _)| name.as_str())
            .collect::<std::collections::HashSet<_>>();
        if wanted.len() < 2 {
            return Ok(Vec::new());
        }

        let directory = tempfile::tempdir()?;
        let candidates = directory.path().join("candidates.fna");
        let written = self.extract(&wanted, &candidates)?;
        info!("Comparing {written} contigs for homology");
        if written < 2 {
            return Ok(Vec::new());
        }
        let table = self.triangle(&candidates)?;
        Ok(parse_table(&table, names))
    }

    fn extract(&self, wanted: &std::collections::HashSet<&str>, path: &Path) -> Result<usize> {
        let mut reader = parse_fastx_file(self.assembly)?;
        let mut sink = BufWriter::new(std::fs::File::create(path)?);
        let mut written = 0;
        while let Some(record) = reader.next() {
            let record = record?;
            let name = String::from_utf8_lossy(record.id())
                .split_whitespace()
                .next()
                .unwrap_or_default()
                .to_string();
            if !wanted.contains(name.as_str()) {
                continue;
            }
            written += 1;
            writeln!(sink, ">{name}")?;
            sink.write_all(&record.seq())?;
            writeln!(sink)?;
        }
        sink.flush()?;
        Ok(written)
    }

    fn triangle(&self, candidates: &Path) -> Result<String> {
        let started = std::time::Instant::now();
        let output = Command::new("skani")
            .args(["triangle", "-i", "-E", "--medium", "--min-af", "0", "-t"])
            .arg(self.threads.to_string())
            .arg(candidates)
            .output()?;
        if !output.status.success() {
            bail!(
                "`skani triangle` failed: {}",
                String::from_utf8_lossy(&output.stderr).trim()
            );
        }
        debug!("skani took {:.1}s", started.elapsed().as_secs_f64());
        Ok(String::from_utf8_lossy(&output.stdout).into_owned())
    }
}

/// skani writes its log to stdout ahead of the table, so the header line is where it starts.
pub fn parse_table(table: &str, names: &[String]) -> Vec<Pair> {
    let index = names
        .iter()
        .enumerate()
        .map(|(row, name)| (name.as_str(), row))
        .collect::<HashMap<_, _>>();
    let mut lines = table
        .lines()
        .skip_while(|line| !line.starts_with("Ref_file"));
    let Some(header) = lines.next() else {
        return Vec::new();
    };
    let column = header
        .split('\t')
        .enumerate()
        .map(|(at, name)| (name, at))
        .collect::<HashMap<_, _>>();
    let wanted = [
        "ANI",
        "Align_fraction_ref",
        "Align_fraction_query",
        "Ref_name",
        "Query_name",
    ]
    .map(|name| column.get(name).copied());
    let [ani, af_ref, af_query, ref_name, query_name] = wanted.map(|at| at.unwrap_or(usize::MAX));

    lines
        .filter_map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            let field = |at: usize| fields.get(at).copied().unwrap_or_default();
            let one = *index.get(field(ref_name).split_whitespace().next()?)?;
            let other = *index.get(field(query_name).split_whitespace().next()?)?;
            if one == other {
                return None;
            }
            Some(Pair {
                one,
                other,
                identity: field(ani).parse().ok()?,
                aligned_one: field(af_ref).parse().ok()?,
                aligned_other: field(af_query).parse().ok()?,
            })
        })
        .collect()
}
