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
            "pass\trung\tworth\tcontigs\tbp\tcompleteness\tcontamination\tverdict\tmembers"
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
        quality: Option<Quality>,
        verdict: &str,
        contigs: &[usize],
    ) {
        let (completeness, contamination) = match quality {
            Some(quality) => (
                format!("{:.2}", quality.completeness),
                format!("{:.2}", quality.contamination),
            ),
            None => ("NA".to_string(), "NA".to_string()),
        };
        let members = contigs
            .iter()
            .filter_map(|contig| self.names.get(*contig))
            .map(String::as_str)
            .collect::<Vec<_>>()
            .join(",");
        let mut sink = self.sink.borrow_mut();
        let _ = writeln!(
            sink,
            "{pass}\t{rung}\t{worth:.4}\t{}\t{bp}\t{completeness}\t{contamination}\t{verdict}\t{members}",
            contigs.len()
        );
    }

    pub fn flush(&self) {
        let _ = self.sink.borrow_mut().flush();
    }
}
