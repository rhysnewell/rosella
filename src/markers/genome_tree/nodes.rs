use std::collections::{HashMap, HashSet};

use anyhow::{Context, Result, bail};

use crate::markers::checkm::bare;

const ROOT_NAME: &str = "root";
const DOMAINS: [&str; 2] = ["k__Bacteria", "k__Archaea"];
// lineage_set hard-codes this floor, since other values break how the selected set is found.
const LEAST_GENOMES: u32 = 2;
const EXCLUDED: [&str; 2] = ["TIGR00398", "TIGR00399"];

pub struct Topology {
    pub parents: Vec<Option<usize>>,
    pub names: Vec<String>,
}

fn ends_name(byte: u8) -> bool {
    matches!(byte, b':' | b',' | b'(' | b')' | b';')
}

pub fn parse_newick(text: &str) -> Result<Topology> {
    let bytes = text.trim().as_bytes();
    let mut parents = vec![None];
    let mut names = vec![String::new()];
    let (mut current, mut at) = (0, 0);
    while at < bytes.len() {
        match bytes[at] {
            byte @ (b'(' | b',') => {
                let parent = match byte {
                    b'(' => current,
                    _ => parents[current].context("the tree has a sibling of its root")?,
                };
                parents.push(Some(parent));
                names.push(String::new());
                current = parents.len() - 1;
                at += 1;
            }
            b')' => {
                current = parents[current].context("the tree closes more clades than it opens")?;
                at += 1;
            }
            b';' => break,
            b':' => {
                at += 1;
                while at < bytes.len() && !ends_name(bytes[at]) {
                    at += 1;
                }
            }
            b'\'' => {
                let length = bytes[at + 1..]
                    .iter()
                    .position(|byte| *byte == b'\'')
                    .context("the tree holds an unclosed quote")?;
                names[current] = String::from_utf8_lossy(&bytes[at + 1..at + 1 + length]).into();
                at += length + 2;
            }
            _ => {
                let start = at;
                while at < bytes.len() && !ends_name(bytes[at]) {
                    at += 1;
                }
                names[current] = String::from_utf8_lossy(&bytes[start..at])
                    .trim()
                    .to_string();
            }
        }
    }
    Ok(Topology { parents, names })
}

pub struct Info {
    pub genomes: u32,
    pub taxonomy: String,
    pub removed: HashSet<String>,
    pub groups: Vec<Vec<String>>,
}

fn listed(field: &str) -> impl Iterator<Item = &str> {
    field.split(',').filter(|entry| !entry.is_empty())
}

pub fn read_nodes(text: &str) -> Result<HashMap<String, Info>> {
    let mut nodes = HashMap::new();
    for line in text.lines().skip(1) {
        let fields = line.split('\t').collect::<Vec<_>>();
        let [uid, genomes, taxonomy, removed, groups] = fields[..] else {
            bail!("a node row holds {} fields, not 5", fields.len());
        };
        let info = Info {
            genomes: genomes.parse()?,
            taxonomy: taxonomy.to_string(),
            removed: listed(removed).map(str::to_string).collect(),
            groups: groups
                .split(';')
                .map(|group| listed(group).map(str::to_string).collect::<Vec<_>>())
                .filter(|group| !group.is_empty())
                .collect(),
        };
        nodes.insert(uid.to_string(), info);
    }
    Ok(nodes)
}

pub struct Chosen {
    pub uid: String,
    pub name: String,
    pub groups: Vec<Vec<String>>,
}

struct Recorded<'a> {
    uid: &'a str,
    name: String,
}

pub struct Tree {
    parents: Vec<Option<usize>>,
    uids: Vec<Option<String>>,
    taxa: Vec<String>,
    leaves: HashMap<String, usize>,
    root: usize,
    info: HashMap<String, Info>,
    selected: HashMap<String, String>,
}

impl Tree {
    pub fn new(newick: &str, nodes: &str, selected: &str) -> Result<Self> {
        let Topology { parents, names } = parse_newick(newick)?;
        let (mut uids, mut taxa, mut leaves) = (Vec::new(), Vec::new(), HashMap::new());
        for (at, name) in names.into_iter().enumerate() {
            let mut fields = name.split('|');
            match name.contains('|') {
                true => {
                    uids.push(fields.next().map(str::to_string));
                    taxa.push(fields.next().unwrap_or_default().to_string());
                }
                false => {
                    uids.push(None);
                    taxa.push(String::new());
                    leaves.insert(name, at);
                }
            }
        }
        let selected = selected
            .lines()
            .skip(1)
            .filter_map(|line| line.split_once('\t'))
            .map(|(uid, chosen)| (uid.to_string(), chosen.trim().to_string()))
            .collect();
        Ok(Self {
            parents,
            uids,
            taxa,
            leaves,
            root: 0,
            info: read_nodes(nodes)?,
            selected,
        })
    }

    pub fn leaf(&self, name: &str) -> Option<usize> {
        self.leaves.get(name).copied()
    }

    pub fn uids(&self) -> impl Iterator<Item = &str> {
        self.uids.iter().flatten().map(String::as_str)
    }

    fn ancestors(&self, node: usize) -> impl Iterator<Item = usize> + '_ {
        std::iter::successors(self.parents[node], |at| self.parents[*at])
    }

    fn named(&self, node: usize, info: &Info) -> String {
        let taxonomy = std::iter::once(info.taxonomy.as_str())
            .chain(
                self.ancestors(node)
                    .filter_map(|at| self.info.get(self.uids[at].as_deref()?))
                    .map(|held| held.taxonomy.as_str()),
            )
            .find(|taxonomy| !taxonomy.is_empty())
            .unwrap_or(ROOT_NAME);
        taxonomy.rsplit(';').next().unwrap_or(ROOT_NAME).to_string()
    }

    // CheckM1 records the root's set under the name root when no node clears its bars.
    fn qualifying(&self, start: usize, force_domain: bool) -> Option<(usize, Recorded<'_>)> {
        let mut last = None;
        let mut node = start;
        loop {
            if let Some(uid) = self.uids[node].as_deref() {
                last = Some(uid);
                let info = self.info.get(uid)?;
                let domain = !force_domain || DOMAINS.contains(&self.taxa[node].as_str());
                if info.genomes >= LEAST_GENOMES && domain {
                    let name = self.named(node, info);
                    return Some((node, Recorded { uid, name }));
                }
            }
            match self.parents[node] {
                Some(parent) => node = parent,
                None => break,
            }
        }
        let name = ROOT_NAME.to_string();
        Some((node, Recorded { uid: last?, name }))
    }

    pub fn choose(&self, leaf: Option<usize>, force_domain: bool) -> Option<Chosen> {
        let Some(leaf) = leaf else {
            let (_, recorded) = self.qualifying(self.root, false)?;
            return self.pick(&[recorded], None);
        };
        let first = self
            .ancestors(leaf)
            .find_map(|at| self.uids[at].as_deref())?;
        let mut recorded = Vec::new();
        let mut current = leaf;
        while let Some(parent) = self.parents[current] {
            let (node, held) = self.qualifying(parent, force_domain)?;
            recorded.push(held);
            current = node;
        }
        self.pick(&recorded, Some(first))
    }

    fn pick(&self, recorded: &[Recorded], first: Option<&str>) -> Option<Chosen> {
        let most = recorded.first()?;
        let mut wanted = self.selected.get(most.uid)?;
        let mut found = None;
        for _ in 0..=self.selected.len() {
            found = recorded.iter().find(|held| held.uid == wanted);
            if found.is_some() {
                break;
            }
            wanted = self.selected.get(wanted)?;
        }
        let found = found?;
        let removed = first
            .and_then(|uid| self.info.get(uid))
            .map(|info| &info.removed);
        let kept = |model: &&String| {
            !EXCLUDED.contains(&model.as_str())
                && removed.is_none_or(|removed| !removed.contains(bare(model)))
        };
        let groups = self
            .info
            .get(found.uid)?
            .groups
            .iter()
            .map(|group| group.iter().filter(kept).cloned().collect::<Vec<_>>())
            .filter(|group| !group.is_empty())
            .collect::<Vec<_>>();
        Some(Chosen {
            uid: found.uid.to_string(),
            name: format!("{} ({})", found.name, found.uid),
            groups,
        })
    }
}
