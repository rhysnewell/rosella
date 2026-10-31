use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::Path;

use anyhow::Result;
use log::debug;

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

pub struct AssemblyGraph {
    pub links: Vec<Link>,
    pub circles: HashSet<String>,
}

// A segment whose only link joins its own end to its own start on one strand is a closed
// circle in the assembler's own graph.
pub fn read_graph<P: AsRef<Path>>(path: P, names: &[String]) -> Result<AssemblyGraph> {
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
    let mut looped = HashSet::new();
    let mut joined = HashSet::new();
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
        if from == to && field(&line, 2) == field(&line, 4) {
            looped.insert(from.to_string());
        } else {
            for end in [from, to] {
                if !joined.contains(end) {
                    joined.insert(end.to_string());
                }
            }
        }
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
    let circles = looped.difference(&joined).cloned().collect::<HashSet<_>>();
    debug!(
        "Assembly graph: {seen} links ({} between kept contigs, {} circles)",
        links.len(),
        circles.len()
    );
    Ok(AssemblyGraph { links, circles })
}
