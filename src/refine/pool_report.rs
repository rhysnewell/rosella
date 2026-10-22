use std::path::Path;

use anyhow::Result;

use crate::quality::Quality;
use crate::report_sink::{Sink, members};

pub struct Row<'a> {
    pub pass: usize,
    pub rung: usize,
    pub worth: f64,
    pub bp: usize,
    pub quality: Quality,
    pub verdict: &'a str,
    pub contigs: &'a [usize],
    pub origins: &'a [(usize, usize)],
}

pub struct PoolReport<'a> {
    names: &'a [String],
    sink: Sink,
}

impl<'a> PoolReport<'a> {
    pub fn create(path: &Path, names: &'a [String]) -> Result<Self> {
        let sink = Sink::create(
            path,
            "pass\trung\tworth\tcontigs\tbp\tcompleteness\tcontamination\tverdict\torigins\tmembers",
        )?;
        Ok(Self { names, sink })
    }

    pub fn row(&self, row: Row<'_>) {
        let Row {
            pass,
            rung,
            worth,
            bp,
            quality,
            verdict,
            contigs,
            origins,
        } = row;
        let completeness = quality.completeness;
        let contamination = quality.contamination;
        let origins = origins
            .iter()
            .map(|(bin, bases)| format!("{bin}:{bases}"))
            .collect::<Vec<_>>()
            .join(",");
        self.sink.line(format_args!(
            "{pass}\t{rung}\t{worth:.4}\t{}\t{bp}\t{completeness:.2}\t{contamination:.2}\t{verdict}\t{origins}\t{}",
            contigs.len(),
            members(self.names, contigs)
        ));
    }

    pub fn flush(&self) {
        self.sink.flush();
    }
}
