use std::cell::RefCell;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

pub struct CombineReport<'a> {
    names: &'a [String],
    lengths: &'a [usize],
    sink: RefCell<BufWriter<File>>,
}

impl<'a> CombineReport<'a> {
    pub fn create(path: &Path, names: &'a [String], lengths: &'a [usize]) -> Result<Self> {
        let mut sink = BufWriter::new(File::create(path)?);
        writeln!(sink, "candidate\tworth\tcontigs\tbp\tverdict\tmembers")?;
        Ok(Self {
            names,
            lengths,
            sink: RefCell::new(sink),
        })
    }

    pub fn row(&self, candidate: usize, worth: f64, verdict: &str, contigs: &[usize]) {
        let bp = contigs
            .iter()
            .filter_map(|contig| self.lengths.get(*contig))
            .sum::<usize>();
        let members = contigs
            .iter()
            .filter_map(|contig| self.names.get(*contig))
            .map(String::as_str)
            .collect::<Vec<_>>()
            .join(",");
        let mut sink = self.sink.borrow_mut();
        let _ = writeln!(
            sink,
            "{candidate}\t{worth:.4}\t{}\t{bp}\t{verdict}\t{members}",
            contigs.len()
        );
    }

    pub fn flush(&self) {
        let _ = self.sink.borrow_mut().flush();
    }
}
