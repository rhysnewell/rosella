use std::collections::HashMap;

use anyhow::Result;
use log::warn;
use std::io::BufRead;

use crate::get_file_reader;

/// Completeness and contamination per bin, read from a CheckM1, CheckM2 or AMBER table.
/// Without one every bin is treated as complete and clean, which is what flight does, so
/// the split decision falls back to internal distance statistics alone.
pub fn read_checkm(path: &str) -> Result<HashMap<String, (f64, f64)>> {
    let reader = get_file_reader(path)?;
    let mut lines = reader.lines();
    let header = match lines.next() {
        Some(header) => header?,
        None => bail!("{} is empty", path),
    };
    let columns = header.trim_end().split('\t').collect::<Vec<_>>();

    let name = ["Bin Id", "Name", "BINID"]
        .iter()
        .find_map(|wanted| columns.iter().position(|column| column == wanted));
    let Some(name) = name else {
        return Err(anyhow!(
            "{} has no Bin Id, Name or BINID column, so its bins cannot be matched up",
            path
        ));
    };

    let completeness = columns.iter().position(|column| *column == "Completeness");
    let contamination = columns.iter().position(|column| *column == "Contamination");
    let precision = columns.iter().position(|column| *column == "precision_bp");
    let recall = columns.iter().position(|column| *column == "recall_bp");

    let mut stats = HashMap::new();
    let mut unreadable = 0;
    for line in lines {
        let line = line?;
        let fields = line.trim_end().split('\t').collect::<Vec<_>>();
        let Some(bin) = fields.get(name) else {
            continue;
        };

        let values = match (completeness, contamination, precision, recall) {
            (Some(complete), Some(contaminated), _, _) => {
                parse_pair(&fields, complete, contaminated)
            }
            // AMBER reports purity and recall as fractions rather than percentages.
            (_, _, Some(precise), Some(recovered)) => parse_pair(&fields, recovered, precise)
                .map(|(recovered, precise)| (recovered * 100.0, (1.0 - precise) * 100.0)),
            _ => {
                return Err(anyhow!(
                    "{} has no Completeness and Contamination columns and no precision_bp and recall_bp pair",
                    path
                ));
            }
        };

        match values {
            Some(values) => {
                stats.insert(bin.to_string(), values);
            }
            None => unreadable += 1,
        }
    }

    if unreadable > 0 {
        warn!("{} rows of {} had unreadable numbers", unreadable, path);
    }

    Ok(stats)
}

fn parse_pair(fields: &[&str], first: usize, second: usize) -> Option<(f64, f64)> {
    let first = fields.get(first)?.trim().parse::<f64>().ok()?;
    let second = fields.get(second)?.trim().parse::<f64>().ok()?;
    Some((first, second))
}
