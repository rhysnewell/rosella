use std::{
    collections::HashMap,
    io::{BufRead, BufWriter, Write},
    path::{Path, PathBuf},
    process::Command,
};

use anyhow::Result;
use log::debug;
use rayon::prelude::*;

pub struct HmmerEngine {
    threads: usize,
    shards: usize,
}

/// HMMER threads a single model block over the sequence database and stops scaling well a few
/// threads in. A protein is scored against every model whatever else is in the file, so cutting
/// the file into pieces and running them at once is the same search.
fn shard_proteins(proteins: &Path, directory: &Path, shards: usize) -> Result<Vec<PathBuf>> {
    let paths = (0..shards)
        .map(|shard| directory.join(format!("shard{shard}.faa")))
        .collect::<Vec<_>>();
    let mut sinks = paths
        .iter()
        .map(|path| Ok(BufWriter::new(std::fs::File::create(path)?)))
        .collect::<Result<Vec<_>>>()?;

    let reader = std::io::BufReader::new(std::fs::File::open(proteins)?);
    let mut at = 0;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            at = (at + 1) % shards;
        }
        writeln!(sinks[at], "{line}")?;
    }
    for sink in &mut sinks {
        sink.flush()?;
    }
    Ok(paths)
}

impl HmmerEngine {
    pub fn new(threads: usize, shards: usize) -> Self {
        Self {
            threads,
            shards: shards.clamp(1, threads.max(1)),
        }
    }

    pub fn check_installed() -> Result<()> {
        match Command::new("hmmsearch").arg("-h").output() {
            Ok(output) if output.status.success() => Ok(()),
            Ok(output) => bail!(
                "`hmmsearch -h` failed: {}",
                String::from_utf8_lossy(&output.stderr).trim()
            ),
            Err(_) => bail!(
                "hmmsearch is not on PATH. rosella judges bins on single copy markers through \
                 HMMER, so install the hmmer package or pass --gene-database instead."
            ),
        }
    }

    /// Best scoring model per protein. A gene that trips two models is one gene, so counting
    /// it under both would inflate presence and duplication at once.
    pub fn search(
        &self,
        hmm: &Path,
        proteins: &Path,
        directory: &Path,
    ) -> Result<HashMap<String, (String, f64)>> {
        let tables = self.run(hmm, proteins, directory, None)?;
        Ok(best_hits(&tables.concat()))
    }

    pub fn search_domains(&self, hmm: &Path, proteins: &Path, directory: &Path) -> Result<String> {
        Ok(self
            .run(hmm, proteins, directory, Some(crate::markers::fragments::DOMAIN_FLOOR))?
            .concat())
    }

    fn run(
        &self,
        hmm: &Path,
        proteins: &Path,
        directory: &Path,
        domains: Option<&str>,
    ) -> Result<Vec<String>> {
        let started = std::time::Instant::now();
        let pieces = match self.shards {
            1 => vec![proteins.to_path_buf()],
            shards => shard_proteins(proteins, directory, shards)?,
        };
        let cpus = (self.threads / self.shards).max(1);

        let tables = pieces
            .par_iter()
            .enumerate()
            .map(|(shard, piece)| {
                let table = directory.join(format!("hits{shard}.tbl"));
                let mut command = Command::new("hmmsearch");
                match domains {
                    Some(floor) => command.args(["--domT", floor, "-T", floor, "--noali", "--cpu"]),
                    None => command.args(["--cut_ga", "--noali", "--cpu"]),
                };
                let output = command
                    .arg(cpus.to_string())
                    .arg("-o")
                    .arg(directory.join(format!("hmmsearch{shard}.log")))
                    .arg(match domains {
                        Some(_) => "--domtblout",
                        None => "--tblout",
                    })
                    .arg(&table)
                    .arg(hmm)
                    .arg(piece)
                    .output()?;
                if !output.status.success() {
                    bail!(
                        "`hmmsearch` failed: {}",
                        String::from_utf8_lossy(&output.stderr).trim()
                    );
                }
                Ok(std::fs::read_to_string(table)?)
            })
            .collect::<Result<Vec<_>>>()?;

        debug!(
            "hmmsearch took {:.1}s over {} shards",
            started.elapsed().as_secs_f64(),
            self.shards
        );
        Ok(tables)
    }
}

pub fn best_hits(table: &str) -> HashMap<String, (String, f64)> {
    let mut best: HashMap<String, (String, f64)> = HashMap::new();
    for line in table.lines().filter(|line| !line.starts_with('#')) {
        let fields = line.split_whitespace().collect::<Vec<_>>();
        let (Some(protein), Some(model), Some(score)) =
            (fields.first(), fields.get(2), fields.get(5))
        else {
            continue;
        };
        let Ok(score) = score.parse::<f64>() else {
            continue;
        };
        match best.get(*protein) {
            Some((_, held)) if *held >= score => {}
            _ => {
                best.insert((*protein).to_string(), ((*model).to_string(), score));
            }
        }
    }
    best
}
