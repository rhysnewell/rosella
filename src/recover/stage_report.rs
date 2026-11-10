use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::Result;
use log::warn;

use crate::clustering::clusterer::Partitioning;

pub struct StageReport {
    directory: PathBuf,
    passes: AtomicUsize,
}

impl StageReport {
    pub fn create(directory: &Path) -> Result<Self> {
        std::fs::create_dir_all(directory)?;
        Ok(Self {
            directory: directory.to_path_buf(),
            passes: AtomicUsize::new(0),
        })
    }

    pub fn next_pass(&self) -> usize {
        self.passes.fetch_add(1, Ordering::Relaxed)
    }

    pub fn partitioning(
        &self,
        name: &str,
        held: &Partitioning,
        contigs: &[usize],
        names: &[String],
    ) {
        let rows = held
            .cluster_map
            .iter()
            .flat_map(|(bin, members)| members.iter().map(move |at| (*bin, contigs[*at])));
        self.write(name, rows, names);
    }

    pub fn bins(&self, name: &str, bins: &BTreeMap<usize, Vec<usize>>, names: &[String]) {
        let rows = bins
            .iter()
            .flat_map(|(bin, members)| members.iter().map(move |contig| (*bin, *contig)));
        self.write(name, rows, names);
    }

    fn write(&self, name: &str, rows: impl Iterator<Item = (usize, usize)>, names: &[String]) {
        let lines = rows.map(|(bin, contig)| format!("{}\t{bin}", names[contig]));
        self.table(name, "contig\tbin", lines);
    }

    pub fn table(&self, name: &str, header: &str, rows: impl Iterator<Item = String>) {
        let path = self.directory.join(format!("{name}.tsv"));
        let written = File::create(&path)
            .map_err(anyhow::Error::from)
            .and_then(|file| {
                let mut out = BufWriter::new(file);
                writeln!(out, "{header}")?;
                for row in rows {
                    writeln!(out, "{row}")?;
                }
                out.flush()?;
                Ok(())
            });
        if let Err(error) = written {
            warn!("Stage report {}: {error}", path.display());
        }
    }
}
