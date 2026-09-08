use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use anyhow::Result;
use log::info;

pub const STAGES_FILE: &str = "stages.tsv";

struct Row {
    stage: &'static str,
    bins: usize,
    binned: usize,
    binned_bp: usize,
    unbinned: usize,
    unbinned_bp: usize,
}

/// Every stage logs its own delta in its own shape, so nothing says where the contigs are
/// after each one. This is the same numbers for every stage, and it needs no truth to read.
#[derive(Default)]
pub struct Census {
    rows: Vec<Row>,
}

impl Census {
    pub fn record<B, C, U>(&mut self, stage: &'static str, bins: B, unbinned: U, lengths: &[usize])
    where
        B: IntoIterator<Item = C>,
        C: IntoIterator<Item = usize>,
        U: IntoIterator<Item = usize>,
    {
        let mut row = Row {
            stage,
            bins: 0,
            binned: 0,
            binned_bp: 0,
            unbinned: 0,
            unbinned_bp: 0,
        };
        for contigs in bins {
            row.bins += 1;
            for contig in contigs {
                row.binned += 1;
                row.binned_bp += lengths[contig];
            }
        }
        for contig in unbinned {
            row.unbinned += 1;
            row.unbinned_bp += lengths[contig];
        }
        self.rows.push(row);
    }

    pub fn write(&self, path: impl AsRef<Path>) -> Result<()> {
        let mut out = BufWriter::new(File::create(path)?);
        writeln!(
            out,
            "stage\tbins\tbinned_contigs\tbinned_bp\tunbinned_contigs\tunbinned_bp"
        )?;
        for row in &self.rows {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}",
                row.stage, row.bins, row.binned, row.binned_bp, row.unbinned, row.unbinned_bp
            )?;
            info!(
                "{} {} bins, {} contigs {} bp binned, {} contigs {} bp unbinned",
                row.stage, row.bins, row.binned, row.binned_bp, row.unbinned, row.unbinned_bp
            );
        }
        out.flush()?;
        Ok(())
    }
}
