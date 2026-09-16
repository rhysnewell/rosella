use std::{
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    process::Command,
};

use anyhow::Result;
use log::{debug, warn};
use rayon::prelude::*;

const SHARD_THREADS: usize = 4;

/// `--cpu n` runs n workers plus a master, so no shard can cost less than this.
const THREADS_PER_SHARD: usize = 2;

const PROTEIN_STEM: &str = "shard";

pub struct HmmerEngine {
    shards: usize,
    cpus: usize,
}

/// HMMER threads a single model block over the sequence database and stops scaling well a few
/// threads in. A protein is scored against every model whatever else is in the file, so cutting
/// the file into pieces and running them at once is the same search.
pub struct Shards {
    sinks: Vec<BufWriter<std::fs::File>>,
    paths: Vec<PathBuf>,
    at: usize,
}

impl Shards {
    pub fn write(&mut self, id: usize, protein: &str) -> Result<()> {
        self.at = (self.at + 1) % self.sinks.len();
        writeln!(self.sinks[self.at], ">{id}\n{protein}")?;
        Ok(())
    }

    pub fn finish(mut self) -> Result<Vec<PathBuf>> {
        for sink in &mut self.sinks {
            sink.flush()?;
        }
        Ok(self.paths)
    }
}

impl HmmerEngine {
    /// The budget is split once here rather than derived from the thread count twice. Four
    /// threads a shard beat both a single wide search and a shard per pair of threads.
    pub fn new(threads: usize, requested: Option<usize>) -> Self {
        let ceiling = (threads / THREADS_PER_SHARD).max(1);
        let shards = requested
            .unwrap_or((threads / SHARD_THREADS).max(1))
            .min(ceiling);
        if let Some(asked) = requested.filter(|asked| *asked > ceiling) {
            warn!("{threads} threads leave room for {ceiling} hmmsearch shards, not {asked}.");
        }
        Self {
            shards,
            cpus: (threads / shards).saturating_sub(1).max(1),
        }
    }

    pub fn protein_shards(&self, directory: &Path) -> Result<Shards> {
        self.open(directory, PROTEIN_STEM)
    }

    fn open(&self, directory: &Path, stem: &str) -> Result<Shards> {
        let paths = (0..self.shards)
            .map(|shard| directory.join(format!("{stem}{shard}.faa")))
            .collect::<Vec<_>>();
        let sinks = paths
            .iter()
            .map(|path| Ok(BufWriter::new(std::fs::File::create(path)?)))
            .collect::<Result<Vec<_>>>()?;
        Ok(Shards {
            sinks,
            paths,
            at: 0,
        })
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
                 HMMER, so install the hmmer package."
            ),
        }
    }

    /// One search serves both readings. The reporting floor sits under the lowest score either
    /// can accept, so the table is a superset of what a gathering-cutoff run would report.
    pub fn search(
        &self,
        hmm: &Path,
        pieces: &[PathBuf],
        directory: &Path,
        floor: &str,
    ) -> Result<String> {
        let started = std::time::Instant::now();
        let progress =
            crate::progress::counted(crate::progress::Stage::SearchingModels, pieces.len() as u64);
        let tables = pieces
            .par_iter()
            .enumerate()
            .map(|(shard, piece)| {
                let table = directory.join(format!("hits{shard}.tbl"));
                let output = Command::new("hmmsearch")
                    .args(["--domT", floor, "-T", floor, "--noali", "--cpu"])
                    .arg(self.cpus.to_string())
                    .arg("-o")
                    .arg(directory.join(format!("hmmsearch{shard}.log")))
                    .arg("--domtblout")
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
                progress.inc(1);
                Ok(std::fs::read_to_string(table)?)
            })
            .collect::<Result<Vec<_>>>()?;
        progress.finish_and_clear();

        debug!(
            "hmmsearch took {:.1}s over {} shards",
            started.elapsed().as_secs_f64(),
            self.shards
        );
        Ok(tables.concat())
    }
}
