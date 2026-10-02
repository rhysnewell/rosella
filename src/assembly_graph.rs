use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::Result;
use log::info;

use crate::get_file_reader;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Link {
    pub from: usize,
    pub to: usize,
}

fn field(line: &[u8], at: usize) -> Option<&str> {
    line.split(|byte| *byte == b'\t')
        .nth(at)
        .and_then(|found| std::str::from_utf8(found).ok())
}

pub fn read_links<P: AsRef<Path>>(path: P, names: &[String]) -> Result<Vec<Link>> {
    let _timer = crate::timing::scope("assembly_graph");
    let index = names
        .iter()
        .enumerate()
        .map(|(at, name)| (name.as_str(), at))
        .collect::<HashMap<_, _>>();
    let mut reader = get_file_reader(path)?;
    let mut line = Vec::new();
    let mut seen = 0;
    let mut links = Vec::new();
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
        if let (Some(a), Some(b)) = (index.get(from), index.get(to))
            && a != b
        {
            links.push(Link {
                from: *a.min(b),
                to: *a.max(b),
            });
        }
    }
    links.sort_unstable();
    links.dedup();
    info!(
        "Assembly graph: {seen} links, {} between contigs that survived the filter.",
        links.len()
    );
    Ok(links)
}
