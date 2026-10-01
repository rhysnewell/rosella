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

// A uniquely named file beside the target, renamed into place once whole, so a killed run never
// leaves a truncated table where a later run will reuse it and two writers never interleave.
pub fn write_atomically(path: &Path, fill: impl FnOnce(&File) -> Result<()>) -> Result<()> {
    let parent = path.parent().unwrap_or(Path::new("."));
    let mut builder = tempfile::Builder::new();
    // A temporary file is owner only, and a marker cache is shared, so the umask decides instead.
    #[cfg(unix)]
    builder.permissions(std::os::unix::fs::PermissionsExt::from_mode(0o666));
    let pending = builder.tempfile_in(parent)?;
    fill(pending.as_file())?;
    pending.persist(path)?;
    Ok(())
}

pub fn members(names: &[String], contigs: &[usize]) -> String {
    contigs
        .iter()
        .filter_map(|contig| names.get(*contig))
        .map(String::as_str)
        .collect::<Vec<_>>()
        .join(",")
}
