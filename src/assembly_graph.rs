use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::Result;
use log::info;

const READ_BUFFER: usize = 1 << 20;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Link {
    pub from: usize,
    pub to: usize,
    pub trust: f32,
}

fn field(line: &[u8], at: usize) -> Option<&str> {
    line.split(|byte| *byte == b'\t')
        .nth(at)
        .and_then(|found| std::str::from_utf8(found).ok())
}

/// A link leaves one end of a segment, and that end's degree is how many continuations the
/// assembler could not choose between. Counted over every link, including those to contigs the
/// length filter dropped, because the branching is there whether or not we kept the neighbour.
#[derive(Default)]
struct Ends {
    ids: HashMap<String, usize>,
    degree: Vec<u32>,
}

impl Ends {
    fn key(&mut self, name: &str, right: bool) -> usize {
        let id = match self.ids.get(name) {
            Some(id) => *id,
            None => {
                let next = self.ids.len();
                self.ids.insert(name.to_string(), next);
                next
            }
        };
        id * 2 + usize::from(right)
    }

    fn count(&mut self, key: usize) {
        if self.degree.len() <= key {
            self.degree.resize(key + 1, 0);
        }
        self.degree[key] += 1;
    }
}

#[derive(Default)]
struct Parsed {
    seen: usize,
    ends: Ends,
    joined: Vec<(usize, usize, usize, usize)>,
    walked: HashSet<(usize, usize)>,
}

fn walk(line: &[u8], index: &HashMap<&str, usize>, walked: &mut HashSet<(usize, usize)>) {
    let Some(steps) = field(line, 2) else {
        return;
    };
    let mut previous: Option<usize> = None;
    for step in steps.split(',') {
        let at = index.get(step.trim_end_matches(['+', '-'])).copied();
        if let (Some(before), Some(now)) = (previous, at)
            && before != now
        {
            walked.insert((before.min(now), before.max(now)));
        }
        previous = at.or(previous);
    }
}

fn parse(path: impl AsRef<Path>, index: &HashMap<&str, usize>) -> Result<Parsed> {
    let mut reader = BufReader::with_capacity(READ_BUFFER, File::open(path)?);
    let mut line = Vec::new();
    let mut parsed = Parsed::default();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        match line.first() {
            Some(b'P') => walk(&line, index, &mut parsed.walked),
            Some(b'L') => {
                parsed.seen += 1;
                let (Some(from), Some(out), Some(to), Some(into)) = (
                    field(&line, 1),
                    field(&line, 2),
                    field(&line, 3),
                    field(&line, 4),
                ) else {
                    continue;
                };
                let leaves = parsed.ends.key(from, out == "+");
                let enters = parsed.ends.key(to, into == "-");
                parsed.ends.count(leaves);
                parsed.ends.count(enters);
                let (Some(a), Some(b)) = (index.get(from), index.get(to)) else {
                    continue;
                };
                if a != b {
                    parsed.joined.push((*a.min(b), *a.max(b), leaves, enters));
                }
            }
            _ => continue,
        }
    }
    Ok(parsed)
}

/// Trust only ever discounts. A link branchier than the assembly's median loses weight and the
/// rest keep the weight the grid settled, because a global raise above it is already refuted.
/// Deriving the pivot from the assembly is what keeps a hub cut off the command line.
fn trusted(
    joined: Vec<(usize, usize, usize, usize)>,
    degree: &[u32],
    walked: &HashSet<(usize, usize)>,
) -> Vec<Link> {
    let mut branching: HashMap<(usize, usize), u32> = HashMap::new();
    for (a, b, leaves, enters) in joined {
        let most = degree[leaves].max(degree[enters]);
        branching
            .entry((a, b))
            .and_modify(|held| *held = (*held).min(most))
            .or_insert(most);
    }
    let mut spread = branching.values().copied().collect::<Vec<_>>();
    if spread.is_empty() {
        return Vec::new();
    }
    spread.sort_unstable();
    let pivot = f32::from(u16::try_from(spread[spread.len() / 2]).unwrap_or(u16::MAX));
    let mut links = branching
        .into_iter()
        .map(|((from, to), most)| {
            let trust = (pivot / f32::from(u16::try_from(most).unwrap_or(u16::MAX))).min(1.0);
            let trust = match walked.contains(&(from, to)) {
                true => 1.0,
                false => trust,
            };
            Link { from, to, trust }
        })
        .collect::<Vec<_>>();
    links.sort_unstable_by_key(|link| (link.from, link.to));
    links
}

pub fn read_links<P: AsRef<Path>>(path: P, names: &[String]) -> Result<Vec<Link>> {
    let _timer = crate::timing::scope("assembly_graph");
    let index = names
        .iter()
        .enumerate()
        .map(|(at, name)| (name.as_str(), at))
        .collect::<HashMap<_, _>>();
    let parsed = parse(path, &index)?;
    let seen = parsed.seen;
    let walked = parsed.walked.len();
    let links = trusted(parsed.joined, &parsed.ends.degree, &parsed.walked);
    let confirmed = links.iter().filter(|link| link.trust >= 1.0).count();
    info!(
        "Assembly graph: {seen} links, {} between contigs that survived the filter, {confirmed} \
         at or above the median branching, {walked} confirmed by a contig path.",
        links.len()
    );
    Ok(links)
}
