use std::path::Path;

use anyhow::Result;

use crate::recover::origin::{Given, row_index};

// Contig groups read from a CAMI binning file, offered to the pool beside its own proposals
// so the bar can be asked whether it would take the right grouping if it were handed one.
pub fn read_groups(path: &str, names: &[String]) -> Result<Vec<Vec<usize>>> {
    let mut groups = Given::table(Path::new(path))?
        .bins_at(&row_index(names), path)
        .into_values()
        .collect::<Vec<_>>();
    groups.sort_unstable();
    Ok(groups)
}
