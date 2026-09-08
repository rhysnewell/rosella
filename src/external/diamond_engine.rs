use std::path::Path;
use std::process::Command;

use anyhow::{Result, bail};

const QUERY_COVER: &str = "80";
const SUBJECT_COVER: &str = "80";
const PERCENT_ID: &str = "30";
const EVALUE: &str = "1e-05";

pub struct DiamondEngine {
    database: std::path::PathBuf,
    threads: usize,
}

impl DiamondEngine {
    pub fn new(database: &Path, threads: usize) -> Result<Self> {
        match Command::new("diamond").arg("--version").output() {
            Ok(output) if output.status.success() => {}
            Ok(output) => bail!(
                "`diamond --version` failed: {}",
                String::from_utf8_lossy(&output.stderr).trim()
            ),
            Err(_) => bail!(
                "diamond is not on PATH. rosella reads gene families through it, so install \
                 the diamond package or drop --checkm2"
            ),
        }
        if !database.is_file() {
            bail!("{} is not a file", database.display());
        }
        Ok(Self {
            database: database.to_path_buf(),
            threads: threads.max(1),
        })
    }

    /// One hit per protein, at the trained model's own thresholds. Anything looser changes the
    /// gene counts the boosters were fitted on.
    pub fn best_hits(&self, proteins: &Path, into: &Path) -> Result<()> {
        let workspace = tempfile::tempdir()?;
        let output = Command::new("diamond")
            .args(["blastp", "--outfmt", "6", "qseqid", "sseqid"])
            .arg("--query")
            .arg(proteins)
            .arg("-o")
            .arg(into)
            .arg("--db")
            .arg(&self.database)
            .args(["--threads", &self.threads.to_string()])
            .args(["--max-target-seqs", "1"])
            .args(["--query-cover", QUERY_COVER])
            .args(["--subject-cover", SUBJECT_COVER])
            .args(["--id", PERCENT_ID])
            .args(["--evalue", EVALUE])
            .arg("--tmpdir")
            .arg(workspace.path())
            .arg("--quiet")
            .output()?;
        if !output.status.success() {
            bail!(
                "diamond blastp failed: {}",
                String::from_utf8_lossy(&output.stderr).trim()
            );
        }
        Ok(())
    }
}
