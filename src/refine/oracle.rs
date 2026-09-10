use std::collections::HashMap;

use anyhow::Result;
use std::io::BufRead;

use crate::get_file_reader;

/// Contig groups read from a CAMI binning file, offered to the pool beside its own proposals
/// so the bar can be asked whether it would take the right grouping if it were handed one.
pub fn read_groups(path: &str, names: &[String]) -> Result<Vec<Vec<usize>>> {
    let index = names
        .iter()
        .enumerate()
        .map(|(position, name)| (name.as_str(), position))
        .collect::<HashMap<_, _>>();

    let mut groups: HashMap<String, Vec<usize>> = HashMap::new();
    let mut unknown = 0;
    for line in get_file_reader(path)?.lines() {
        let line = line?;
        if line.starts_with(['@', '#']) || line.trim().is_empty() {
            continue;
        }
        let mut fields = line.trim_end().split('\t');
        let (Some(contig), Some(group)) = (fields.next(), fields.next()) else {
            continue;
        };
        if contig == "SEQUENCEID" {
            continue;
        }
        match index.get(contig) {
            Some(position) => groups.entry(group.to_string()).or_default().push(*position),
            None => unknown += 1,
        }
    }

    if unknown > 0 {
        log::warn!("{unknown} contigs of {path} are not in this assembly");
    }
    let mut groups = groups
        .into_values()
        .map(|mut contigs| {
            contigs.sort_unstable();
            contigs
        })
        .collect::<Vec<_>>();
    groups.sort_unstable();
    Ok(groups)
}
