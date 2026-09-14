use std::cell::RefCell;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::quality::Quality;

pub struct PoolReport<'a> {
    names: &'a [String],
    sink: RefCell<BufWriter<File>>,
}

impl<'a> PoolReport<'a> {
    pub fn create(path: &Path, names: &'a [String]) -> Result<Self> {
        let mut sink = BufWriter::new(File::create(path)?);
        writeln!(
            sink,
            "pass\trung\tworth\tcontigs\tbp\tcompleteness\tcontamination\tverdict\torigins\tmembers"
        )?;
        Ok(Self {
            names,
            sink: RefCell::new(sink),
        })
    }

    pub fn row(
        &self,
        pass: usize,
        rung: usize,
        worth: f64,
        bp: usize,
        quality: Quality,
        verdict: &str,
        contigs: &[usize],
        origins: &[(usize, usize)],
    ) {
        let completeness = quality.completeness;
        let contamination = quality.contamination;
        let members = contigs
            .iter()
            .filter_map(|contig| self.names.get(*contig))
            .map(String::as_str)
            .collect::<Vec<_>>()
            .join(",");
        let origins = origins
            .iter()
            .map(|(bin, bases)| format!("{bin}:{bases}"))
            .collect::<Vec<_>>()
            .join(",");
        let mut sink = self.sink.borrow_mut();
        let _ = writeln!(
            sink,
            "{pass}\t{rung}\t{worth:.4}\t{}\t{bp}\t{completeness:.2}\t{contamination:.2}\t{verdict}\t{origins}\t{members}",
            contigs.len()
        );
    }

    pub fn flush(&self) {
        let _ = self.sink.borrow_mut().flush();
    }
}
