use std::collections::HashMap;

use crate::kmers::sketch::ContigSketches;
use crate::refine::duplication::{DEFAULT_MIN_HASHES, MAX_SPREAD};

pub const DEFAULT_APART: f64 = 0.5;
pub const DEFAULT_TOGETHER: f64 = 0.9;
pub const LINK_SCOPE_NAMES: [&str; 3] = ["both", "apart", "together"];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LinkScope {
    Both,
    Apart,
    Together,
}

impl LinkScope {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "both" => Some(Self::Both),
            "apart" => Some(Self::Apart),
            "together" => Some(Self::Together),
            _ => None,
        }
    }

    fn holds_apart(self) -> bool {
        self != Self::Together
    }

    fn holds_together(self) -> bool {
        self != Self::Apart
    }
}

#[derive(Debug, Clone, Copy)]
pub struct LinkSettings {
    pub min_hashes: usize,
    pub apart: f64,
    pub together: f64,
    pub scope: LinkScope,
}

impl Default for LinkSettings {
    fn default() -> Self {
        Self {
            min_hashes: DEFAULT_MIN_HASHES,
            apart: DEFAULT_APART,
            together: DEFAULT_TOGETHER,
            scope: LinkScope::Both,
        }
    }
}

#[derive(Debug, Default)]
pub struct Links {
    pub apart: Vec<(usize, usize)>,
    pub together: Vec<(usize, usize)>,
}

impl Links {
    pub fn components(&self, contigs: usize) -> Vec<usize> {
        let mut parent = (0..contigs).collect::<Vec<_>>();
        for (one, other) in &self.together {
            let (left, right) = (find(&mut parent, *one), find(&mut parent, *other));
            if left != right {
                parent[left] = right;
            }
        }
        (0..contigs).map(|at| find(&mut parent, at)).collect()
    }
}

fn find(parent: &mut [usize], mut at: usize) -> usize {
    while parent[at] != at {
        parent[at] = parent[parent[at]];
        at = parent[at];
    }
    at
}

/// The assembler already merged what it could, so near identical over both lengths is one locus in
/// two organisms, while held whole inside one longer contig and nowhere else is one genome twice.
pub fn links(sketches: &ContigSketches, settings: LinkSettings) -> Links {
    let mut owners: HashMap<u64, Vec<u32>> = HashMap::new();
    for contig in 0..sketches.len() {
        for hash in sketches.hashes(contig) {
            owners.entry(*hash).or_default().push(contig as u32);
        }
    }

    let mut shared: HashMap<(u32, u32), u32> = HashMap::new();
    for holders in owners.values() {
        if holders.len() < 2 || holders.len() > MAX_SPREAD {
            continue;
        }
        for left in 0..holders.len() {
            for right in left + 1..holders.len() {
                *shared.entry((holders[left], holders[right])).or_default() += 1;
            }
        }
    }

    let mut partners = vec![0u32; sketches.len()];
    let mut contained = Vec::new();
    let mut links = Links::default();
    for ((one, other), count) in &shared {
        let (one, other) = (*one as usize, *other as usize);
        let (held_one, held_other) = (sketches.hashes(one).len(), sketches.hashes(other).len());
        if held_one < settings.min_hashes || held_other < settings.min_hashes {
            continue;
        }
        let reach_one = f64::from(*count) / held_one as f64;
        let reach_other = f64::from(*count) / held_other as f64;
        if reach_one.min(reach_other) >= settings.apart {
            if settings.scope.holds_apart() {
                links.apart.push((one, other));
            }
            continue;
        }
        if reach_one.max(reach_other) < settings.together || !settings.scope.holds_together() {
            continue;
        }
        let inside = if reach_one > reach_other { one } else { other };
        partners[inside] += 1;
        contained.push((one, other, inside));
    }

    links.together = contained
        .into_iter()
        .filter(|(_, _, inside)| partners[*inside] == 1)
        .map(|(one, other, _)| (one, other))
        .collect();
    links.apart.sort_unstable();
    links.together.sort_unstable();
    links
}
