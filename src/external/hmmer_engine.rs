use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::{Result, anyhow};
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

pub fn align(profile: &str, proteins: &[&str]) -> Result<String> {
    super::hmmer_library::align(profile, proteins)
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
        "HMMER: {} models, {slices} slices, {:.1}s",
        models.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(table)
}
