use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use log::warn;

use crate::recover::recover_engine::RecoverEngine;
use crate::tables::Tables;

impl RecoverEngine {
    pub(super) fn write_points(&self) {
        let Some(path) = &self.reports.points_report else {
            return;
        };
        if self.tables.learned.is_none() {
            warn!("Points report: no learned space");
            return;
        }
        if let Err(error) = write(path, &self.tables) {
            warn!("Points report: {error} ({})", path.display());
        }
    }
}

fn write(path: &Path, tables: &Tables) -> Result<()> {
    let points = &tables.tnf.kmer_table;
    let mut out = BufWriter::new(File::create(path)?);
    let columns = (1..=points.ncols())
        .map(|at| format!("d{at}"))
        .collect::<Vec<_>>();
    writeln!(out, "contig\tlength\t{}", columns.join("\t"))?;
    for (row, point) in points.rows().into_iter().enumerate() {
        write!(
            out,
            "{}\t{}",
            tables.coverage.contig_names[row], tables.coverage.contig_lengths[row]
        )?;
        for value in point {
            write!(out, "\t{value:.6}")?;
        }
        writeln!(out)?;
    }
    out.flush()?;
    Ok(())
}
