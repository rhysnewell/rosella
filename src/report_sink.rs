use std::fs::File;
use std::path::Path;

use anyhow::Result;

// A uniquely named file beside the target, renamed into place once whole, so a killed run never
// leaves a truncated table where a later run will reuse it and two writers never interleave.
pub fn write_atomically(path: &Path, fill: impl FnOnce(&File) -> Result<()>) -> Result<()> {
    let pending = pending_beside(path)?;
    fill(pending.as_file())?;
    pending.persist(path)?;
    Ok(())
}

// Dropped unpersisted on any failure, so a failed step leaves nothing behind.
pub fn pending_beside(path: &Path) -> Result<tempfile::NamedTempFile> {
    let parent = path.parent().unwrap_or(Path::new("."));
    let mut builder = tempfile::Builder::new();
    // A temporary file is owner only, and a marker cache is shared, so the umask decides instead.
    #[cfg(unix)]
    builder.permissions(std::os::unix::fs::PermissionsExt::from_mode(0o666));
    Ok(builder.tempfile_in(parent)?)
}
