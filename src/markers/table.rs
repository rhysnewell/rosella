use std::collections::{HashMap, HashSet};

use crate::markers::panels::Panel;
use crate::markers::{checkm, sets};

const TABLE: &str = include_str!("../../data/gtdb_markers.tsv");
const SET_TABLE: &str = include_str!("../../data/marker_sets.tsv");

pub struct MarkerSet {
    ids: HashMap<String, u16>,
    names: Vec<String>,
    accessions: HashSet<String>,
    pub(crate) sets: sets::Sets,
    pub(crate) checkm: checkm::Panel,
    cutoffs: Vec<Option<f32>>,
}

impl MarkerSet {
    pub fn embedded() -> Self {
        Self::joined(None)
    }

    pub fn joined(panel: Option<&Panel>) -> Self {
        let mut set = Self::parse_with(TABLE, panel).with_bounds(SET_TABLE);
        set.checkm = checkm::Panel::embedded(|model| set.id(model));
        set
    }

    pub fn parse(table: &str) -> Self {
        Self::parse_with(table, None)
    }

    pub fn parse_with(table: &str, panel: Option<&Panel>) -> Self {
        let mut lines = table.lines();
        let header = lines
            .next()
            .unwrap_or_default()
            .split('\t')
            .collect::<Vec<_>>();
        let column = |name: &str| header.iter().position(|field| *field == name);
        let Some(name_at) = column("model_name") else {
            return Self {
                ids: HashMap::new(),
                names: Vec::new(),
                accessions: HashSet::new(),
                sets: sets::Sets::default(),
                checkm: checkm::Panel::default(),
                cutoffs: Vec::new(),
            };
        };
        let (set_at, accession_at) = (column("sets"), column("accession"));
        let mut group_names = set_names(table);
        let rate_at = group_names
            .iter()
            .map(|group| column(&format!("ubiquity_{group}")))
            .collect::<Vec<_>>();

        let mut ids = HashMap::new();
        let mut names = Vec::new();
        let mut accessions = HashSet::new();
        let mut member_of = vec![Vec::new(); group_names.len()];
        let mut rates = vec![Vec::new(); group_names.len()];
        for line in lines {
            let fields = line.split('\t').collect::<Vec<_>>();
            let Some(name) = fields.get(name_at) else {
                continue;
            };
            if ids.contains_key(*name) {
                continue;
            }
            ids.insert((*name).to_string(), names.len() as u16);
            names.push((*name).to_string());
            if let Some(accession) = accession_at.and_then(|at| fields.get(at)) {
                accessions.insert(checkm::bare(accession).to_string());
            }
            for (group, held) in group_names.iter().enumerate() {
                let member = set_at
                    .and_then(|at| fields.get(at))
                    .is_some_and(|listed| listed.split(',').any(|entry| entry == held));
                member_of[group].push(member);
                let rate = rate_at[group]
                    .and_then(|at| fields.get(at))
                    .and_then(|field| field.parse::<f64>().ok())
                    .unwrap_or(f64::from(u8::from(member)));
                rates[group].push(rate);
            }
        }
        let mut cutoffs = vec![None; names.len()];
        let prokaryotic = group_names.len();
        for lineage in panel.iter().flat_map(|panel| &panel.lineages) {
            let group = group_names.len();
            group_names.push(lineage.name.to_string());
            member_of.push(vec![false; names.len()]);
            rates.push(vec![0.0; names.len()]);
            for (model, cutoff) in &lineage.models {
                if ids.contains_key(model) {
                    continue;
                }
                ids.insert(model.clone(), names.len() as u16);
                names.push(model.clone());
                cutoffs.push(Some(*cutoff as f32));
                for (held, (member, rate)) in member_of.iter_mut().zip(&mut rates).enumerate() {
                    member.push(held == group);
                    rate.push(f64::from(u8::from(held == group)));
                }
            }
        }
        let eukaryotic = (0..group_names.len())
            .map(|group| group >= prokaryotic)
            .collect();
        Self {
            ids,
            names,
            accessions,
            sets: sets::Sets::new(group_names, member_of, &rates).with_eukaryotic(eukaryotic),
            checkm: checkm::Panel::default(),
            cutoffs,
        }
    }

    pub fn holds_family(&self, family: &str) -> bool {
        self.accessions.contains(family)
    }

    // A eukaryotic call short of the whole gene's cutoff is one exon of it, never a second copy.
    pub fn piece(&self, marker: u16, score: f32) -> bool {
        self.cutoff(marker).is_some_and(|cutoff| score < cutoff)
    }

    pub fn has_panel(&self) -> bool {
        self.cutoffs.last().is_some_and(Option::is_some)
    }

    pub fn eukaryotic(&self, marker: u16) -> bool {
        self.cutoff(marker).is_some()
    }

    fn cutoff(&self, marker: u16) -> Option<f32> {
        self.cutoffs.get(marker as usize).copied().flatten()
    }

    pub fn with_bounds(mut self, table: &str) -> Self {
        let mut bounds = vec![0.0; self.sets.len()];
        let mut lines = table.lines();
        let header = lines
            .next()
            .unwrap_or_default()
            .split('\t')
            .collect::<Vec<_>>();
        let column = |name: &str| header.iter().position(|field| *field == name);
        let (Some(set_at), Some(max_at)) = (column("set"), column("max_genome_bp")) else {
            return self;
        };
        for line in lines {
            let fields = line.split('\t').collect::<Vec<_>>();
            let Some(found) = fields
                .get(set_at)
                .and_then(|name| (0..self.sets.len()).find(|set| self.sets.name(*set) == *name))
            else {
                continue;
            };
            if let Some(bp) = fields.get(max_at).and_then(|field| field.parse().ok()) {
                bounds[found] = bp;
            }
        }
        self.sets = self.sets.with_bounds(bounds);
        self
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    pub fn id(&self, model: &str) -> Option<u16> {
        self.ids.get(model).copied()
    }

    pub fn name(&self, marker: u16) -> &str {
        self.names
            .get(marker as usize)
            .map(String::as_str)
            .unwrap_or_default()
    }
}

fn set_names(table: &str) -> Vec<String> {
    let mut lines = table.lines();
    let header = lines
        .next()
        .unwrap_or_default()
        .split('\t')
        .collect::<Vec<_>>();
    let Some(at) = header.iter().position(|field| *field == "sets") else {
        return Vec::new();
    };
    let mut found = Vec::new();
    for line in lines {
        let Some(listed) = line.split('\t').nth(at) else {
            continue;
        };
        for name in listed.split(',').filter(|name| !name.is_empty()) {
            if !found.iter().any(|held| held == name) {
                found.push(name.to_string());
            }
        }
    }
    found
}
