use std::cell::RefCell;
use std::fmt::Arguments;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

pub struct Sink {
    writer: RefCell<BufWriter<File>>,
}

impl Sink {
    pub fn create(path: &Path, header: &str) -> Result<Self> {
        let mut writer = BufWriter::new(File::create(path)?);
        writeln!(writer, "{header}")?;
        Ok(Self {
            writer: RefCell::new(writer),
        })
    }

    // A report is a side channel, so a row that will not write is dropped rather than failing
    // the run that produced it.
    pub fn line(&self, row: Arguments<'_>) {
        let _ = writeln!(self.writer.borrow_mut(), "{row}");
    }

    pub fn flush(&self) {
        let _ = self.writer.borrow_mut().flush();
    }
}

pub fn members(names: &[String], contigs: &[usize]) -> String {
    contigs
        .iter()
        .filter_map(|contig| names.get(*contig))
        .map(String::as_str)
        .collect::<Vec<_>>()
        .join(",")
}
