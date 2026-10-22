use std::path::Path;

use anyhow::Result;

use crate::report_sink::{Sink, members};

pub struct CombineReport<'a> {
    names: &'a [String],
    lengths: &'a [usize],
    sink: Sink,
}

impl<'a> CombineReport<'a> {
    pub fn create(path: &Path, names: &'a [String], lengths: &'a [usize]) -> Result<Self> {
        let sink = Sink::create(path, "candidate\tworth\tcontigs\tbp\tverdict\tmembers")?;
        Ok(Self {
            names,
            lengths,
            sink,
        })
    }

    pub fn row(&self, candidate: usize, worth: f64, verdict: &str, contigs: &[usize]) {
        let bp = contigs
            .iter()
            .filter_map(|contig| self.lengths.get(*contig))
            .sum::<usize>();
        self.sink.line(format_args!(
            "{candidate}\t{worth:.4}\t{}\t{bp}\t{verdict}\t{}",
            contigs.len(),
            members(self.names, contigs)
        ));
    }

    pub fn flush(&self) {
        self.sink.flush();
    }
}
