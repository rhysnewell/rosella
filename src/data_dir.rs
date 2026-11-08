use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

const VARIABLE: &str = "ROSELLA_DATA_DIR";

// An active conda or pixi environment keeps its own record, so each environment finds its own download.
fn record() -> Option<PathBuf> {
    let base = std::env::var_os("CONDA_PREFIX")
        .map(|prefix| PathBuf::from(prefix).join("etc"))
        .or_else(|| std::env::var_os("XDG_CONFIG_HOME").map(PathBuf::from))
        .or_else(|| std::env::var_os("HOME").map(|home| PathBuf::from(home).join(".config")))?;
    Some(base.join("rosella").join("data_dir"))
}

pub fn remembered() -> Option<PathBuf> {
    if let Some(directory) = std::env::var_os(VARIABLE) {
        return Some(PathBuf::from(directory));
    }
    let text = std::fs::read_to_string(record()?).ok()?;
    let directory = text.trim();
    (!directory.is_empty()).then(|| PathBuf::from(directory))
}

pub fn remember(directory: &Path) -> Result<()> {
    let path = record()
        .with_context(|| format!("no CONDA_PREFIX, XDG_CONFIG_HOME or HOME to remember the data directory in, so set {VARIABLE}"))?;
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let directory = std::fs::canonicalize(directory)?;
    std::fs::write(&path, format!("{}\n", directory.display()))
        .with_context(|| format!("could not write {}", path.display()))
}
