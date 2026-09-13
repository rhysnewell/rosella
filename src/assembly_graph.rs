use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::Result;
use log::info;

const READ_BUFFER: usize = 1 << 20;

fn field(line: &[u8], at: usize) -> Option<&str> {
    line.split(|byte| *byte == b'\t')
        .nth(at)
        .and_then(|found| std::str::from_utf8(found).ok())
}

pub fn read_links<P: AsRef<Path>>(path: P, names: &[String]) -> Result<Vec<(usize, usize)>> {
    let _timer = crate::timing::scope("assembly_graph");
    let index = names
        .iter()
        .enumerate()
        .map(|(at, name)| (name.as_str(), at))
        .collect::<HashMap<_, _>>();
    let mut reader = BufReader::with_capacity(READ_BUFFER, File::open(path)?);
    let mut line = Vec::new();
    let mut found = Vec::new();
    let mut seen = 0usize;
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.first() != Some(&b'L') {
            continue;
        }
        seen += 1;
        let (Some(from), Some(to)) = (field(&line, 1), field(&line, 3)) else {
            continue;
        };
        let (Some(a), Some(b)) = (index.get(from), index.get(to)) else {
            continue;
        };
        if a != b {
            found.push((*a.min(b), *a.max(b)));
        }
    }
    found.sort_unstable();
    found.dedup();
    info!(
        "Assembly graph: {seen} links, {} between contigs that survived the filter.",
        found.len()
    );
    Ok(found)
}
