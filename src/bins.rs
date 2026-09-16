use std::path::{Path, PathBuf};

use anyhow::{Result, bail};

/// `-f` and `-d` are one arg group, so a run naming both means both. The extension is
/// normalised here because a user who writes `-x .fna` means the same thing as `-x fna`.
pub fn discover(
    files: &[String],
    directory: Option<&String>,
    extension: &str,
) -> Result<Vec<PathBuf>> {
    let wanted = extension.trim_start_matches('.');
    let mut found = files.iter().map(PathBuf::from).collect::<Vec<_>>();
    if let Some(directory) = directory {
        for entry in std::fs::read_dir(directory)? {
            let path = entry?.path();
            if path.extension().and_then(|found| found.to_str()) == Some(wanted) {
                found.push(path);
            }
        }
    }
    found.sort();
    found.dedup();
    if found.is_empty() {
        bail!("no bins found. Pass them with --genome-fasta-files or --genome-fasta-directory");
    }
    Ok(found)
}

pub fn stem(path: &Path) -> String {
    path.file_stem()
        .and_then(|stem| stem.to_str())
        .unwrap_or("bin")
        .to_string()
}

/// A second run into a used directory would append to the bin files already there, so a
/// directory holding output refuses the run rather than corrupting it.
pub fn refuse_used(directory: &str) -> Result<()> {
    let path = Path::new(directory);
    if !path.exists() {
        return Ok(());
    }
    for entry in std::fs::read_dir(path)? {
        let entry = entry?;
        if entry.path().extension().and_then(|found| found.to_str())
            == Some(crate::defaults::FASTA_EXTENSION)
        {
            bail!(
                "{directory} already holds .{} files. Remove them or name another output directory",
                crate::defaults::FASTA_EXTENSION
            );
        }
    }
    Ok(())
}
