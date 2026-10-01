use std::collections::HashMap;

use anyhow::Result;
use log::warn;
use std::io::BufRead;

use crate::get_file_reader;

/// Without a table every bin is clean, so the split falls back to distances.
pub fn read_contamination(path: &str) -> Result<HashMap<String, f64>> {
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

    // AMBER reports purity as a fraction rather than contamination as a percentage.
    let position = |wanted: &str| columns.iter().position(|column| *column == wanted);
    let (column, purity) = match (position("Contamination"), position("precision_bp")) {
        (Some(at), _) => (at, false),
        (None, Some(at)) => (at, true),
        (None, None) => bail!("{path} has no Contamination column and no precision_bp column"),
    };

    let mut stats = HashMap::new();
    let mut unreadable = 0;
    for line in lines {
        let line = line?;
        let fields = line.trim_end().split('\t').collect::<Vec<_>>();
        let Some(bin) = fields.get(name) else {
            continue;
        };
        match fields
            .get(column)
            .and_then(|field| field.trim().parse::<f64>().ok())
        {
            Some(value) if purity => {
                stats.insert(bin.to_string(), (1.0 - value) * 100.0);
            }
            Some(value) => {
                stats.insert(bin.to_string(), value);
            }
            None => unreadable += 1,
        }
    }

    if unreadable > 0 {
        warn!("{} rows of {} had unreadable numbers", unreadable, path);
    }

    Ok(stats)
}
