use std::path::Path;
use std::process::Command;

use anyhow::{Result, bail};

const QUERY_COVER: &str = "80";
const SUBJECT_COVER: &str = "80";
const PERCENT_ID: &str = "30";
const EVALUE: &str = "1e-05";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sensitivity {
    Default,
    Fast,
    Faster,
}

impl Sensitivity {
    pub fn parse(value: &str) -> Option<Self> {
        match value {
            "default" => Some(Self::Default),
            "fast" => Some(Self::Fast),
            "faster" => Some(Self::Faster),
            _ => None,
        }
    }

    pub fn tag(&self) -> &'static str {
        match self {
            Self::Default => "default",
            Self::Fast => "fast",
            Self::Faster => "faster",
        }
    }

    fn flag(&self) -> Option<&'static str> {
        match self {
            Self::Default => None,
            Self::Fast => Some("--fast"),
            Self::Faster => Some("--faster"),
        }
    }
}

pub struct DiamondEngine {
    database: std::path::PathBuf,
    threads: usize,
    sensitivity: Sensitivity,
}

impl DiamondEngine {
    pub fn new(database: &Path, threads: usize, sensitivity: Sensitivity) -> Result<Self> {
        match Command::new("diamond").arg("--version").output() {
            Ok(output) if output.status.success() => {}
            Ok(output) => bail!(
                "`diamond --version` failed: {}",
                String::from_utf8_lossy(&output.stderr).trim()
            ),
            Err(_) => bail!(
                "diamond is not on PATH. rosella reads gene families through it, so install \
                 the diamond package or run without a gene family database"
            ),
        }
        if !database.is_file() {
            bail!("{} is not a file", database.display());
        }
        Ok(Self {
            database: database.to_path_buf(),
            threads: threads.max(1),
            sensitivity,
        })
    }

    /// One hit per protein, at the trained model's own thresholds. Anything looser changes the
    /// gene counts the boosters were fitted on.
    pub fn best_hits(&self, proteins: &Path, into: &Path, workspace: &Path) -> Result<()> {
        let mut command = Command::new("diamond");
        command
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
            .arg(workspace)
            .arg("--quiet");
        if let Some(tier) = self.sensitivity.flag() {
            command.arg(tier);
        }
        let output = command.output()?;
        if !output.status.success() {
            bail!(
                "diamond blastp failed: {}",
                String::from_utf8_lossy(&output.stderr).trim()
            );
        }
        Ok(())
    }
}
