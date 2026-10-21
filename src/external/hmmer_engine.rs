use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    process::Command,
};

use anyhow::{Context, Result, anyhow, bail};
use log::debug;
use rayon::prelude::*;

use super::hmmer_library::{Models, Targets};

// The proteins are searched a slice at a time, so a large assembly never holds them all at once.
const HELD_RESIDUES: usize = 1 << 28;

pub struct ProteinFile {
    sink: BufWriter<File>,
    path: PathBuf,
    filled: bool,
}

impl ProteinFile {
    pub fn create(path: PathBuf) -> Result<Self> {
        Ok(Self {
            sink: BufWriter::new(File::create(&path)?),
            path,
            filled: false,
        })
    }

    pub fn write(&mut self, id: usize, protein: &str) -> Result<()> {
        writeln!(self.sink, ">{id}\n{protein}")?;
        self.filled = true;
        Ok(())
    }

    pub fn finish(mut self) -> Result<Option<PathBuf>> {
        self.sink.flush()?;
        Ok(self.filled.then_some(self.path))
    }
}

pub fn each_protein(path: &Path, mut visit: impl FnMut(usize, &str) -> Result<()>) -> Result<()> {
    let mut lines = BufReader::new(File::open(path)?).lines();
    while let (Some(header), Some(protein)) = (lines.next(), lines.next()) {
        let header = header?;
        let id = header
            .strip_prefix('>')
            .and_then(|id| id.parse::<usize>().ok())
            .ok_or_else(|| anyhow!("{} holds a bad header {header}", path.display()))?;
        visit(id, &protein?)?;
    }
    Ok(())
}

// A run reads and embeds every contig before its first alignment, so a missing hmmalign has to
// stop it at the door rather than minutes in.
pub fn check_aligner() -> Result<()> {
    match Command::new("hmmalign").arg("-h").output() {
        Ok(output) if output.status.success() => Ok(()),
        Ok(output) => bail!(
            "`hmmalign -h` failed: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        ),
        Err(_) => bail!(
            "hmmalign is not on PATH. rosella aligns CheckM's markers for strain heterogeneity \
             through HMMER, so install the hmmer package."
        ),
    }
}

pub fn align(hmm: &Path, sequences: &Path) -> Result<String> {
    let output = Command::new("hmmalign")
        .args(["--outformat", "A2M"])
        .arg(hmm)
        .arg(sequences)
        .output()
        .context("running hmmalign from the hmmer package")?;
    if !output.status.success() {
        bail!(
            "`hmmalign` failed: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        );
    }
    Ok(String::from_utf8(output.stdout)?)
}

pub fn search(hmm: &Path, proteins: &Path, floor: f64) -> Result<String> {
    search_held(hmm, proteins, floor, HELD_RESIDUES)
}

// One search serves both readings. The reporting floor sits under the lowest score either can
// accept, so the table is a superset of what a gathering-cutoff run would report.
pub fn search_held(hmm: &Path, proteins: &Path, floor: f64, held: usize) -> Result<String> {
    let started = std::time::Instant::now();
    let models = Models::read(hmm)?;
    let mut targets = Targets::new()?;
    let progress =
        crate::progress::counted(crate::progress::Stage::SearchingModels, models.len() as u64);
    let mut table = String::new();
    let mut slices = 0;
    let mut slice = |targets: &mut Targets| -> Result<()> {
        slices += 1;
        progress.set_length((slices * models.len()) as u64);
        let tables = crate::pool::install(|| {
            (0..models.len())
                .into_par_iter()
                .map(|model| {
                    let found = models.search(model, targets, floor);
                    progress.inc(1);
                    found
                })
                .collect::<Result<Vec<_>>>()
        })?;
        table.extend(tables);
        targets.clear();
        Ok(())
    };
    each_protein(proteins, |id, protein| {
        targets.add(&models, id, protein)?;
        match targets.residues() >= held {
            true => slice(&mut targets),
            false => Ok(()),
        }
    })?;
    if targets.residues() > 0 {
        slice(&mut targets)?;
    }
    crate::progress::done(&progress);
    debug!(
        "HMMER searched {} models in {slices} slices in {:.1}s",
        models.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(table)
}
