use std::{collections::HashMap, path::Path, process::Command};

use anyhow::Result;
use log::debug;

pub struct HmmerEngine {
    threads: usize,
}

impl HmmerEngine {
    pub fn new(threads: usize) -> Self {
        Self { threads }
    }

    pub fn check_installed() -> Result<()> {
        match Command::new("hmmsearch").arg("-h").output() {
            Ok(output) if output.status.success() => Ok(()),
            Ok(output) => bail!(
                "`hmmsearch -h` failed: {}",
                String::from_utf8_lossy(&output.stderr).trim()
            ),
            Err(_) => bail!(
                "hmmsearch is not on PATH. rosella finds single copy markers through HMMER, so \
                 install the hmmer package or drop --markers"
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
        let started = std::time::Instant::now();
        let table = directory.join("hits.tbl");
        let output = Command::new("hmmsearch")
            .args(["--cut_ga", "--noali", "--cpu"])
            .arg(self.threads.to_string())
            .arg("-o")
            .arg(directory.join("hmmsearch.log"))
            .arg("--tblout")
            .arg(&table)
            .arg(hmm)
            .arg(proteins)
            .output()?;
        if !output.status.success() {
            bail!(
                "`hmmsearch` failed: {}",
                String::from_utf8_lossy(&output.stderr).trim()
            );
        }
        debug!("hmmsearch took {:.1}s", started.elapsed().as_secs_f64());
        Ok(best_hits(&std::fs::read_to_string(table)?))
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
