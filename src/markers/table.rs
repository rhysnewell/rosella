use std::collections::HashMap;
use std::ops::Range;

use crate::markers::panels::Panel;
use crate::markers::{checkm, sets};

const TABLE: &str = include_str!("../../data/gtdb_markers.tsv");
const SET_TABLE: &str = include_str!("../../data/marker_sets.tsv");
const UBIQUITY: &str = include_str!("../../data/checkm_ubiquity.tsv");

pub const FAMILY_PREFIX: &str = "checkm_";

pub struct MarkerSet {
    ids: HashMap<String, u16>,
    names: Vec<String>,
    pub(crate) sets: sets::Sets,
    pub(crate) checkm: checkm::Panel,
    cutoffs: Vec<Option<f32>>,
    families: Vec<Option<u16>>,
    from_checkm: Range<usize>,
}

#[derive(Default)]
struct Draft {
    ids: HashMap<String, u16>,
    names: Vec<String>,
    accessions: HashMap<String, u16>,
    groups: Vec<String>,
    member_of: Vec<Vec<bool>>,
    rates: Vec<Vec<f64>>,
    cutoffs: Vec<Option<f32>>,
    families: Vec<Option<u16>>,
    from_checkm: Range<usize>,
    prokaryotic: Option<usize>,
}

impl Draft {
    fn parse(table: &str) -> Self {
        let mut lines = table.lines();
        let header = lines
            .next()
            .unwrap_or_default()
            .split('\t')
            .collect::<Vec<_>>();
        let column = |name: &str| header.iter().position(|field| *field == name);
        let Some(name_at) = column("model_name") else {
            return Self::default();
        };
        let (set_at, accession_at) = (column("sets"), column("accession"));
        let groups = set_names(table);
        let rate_at = groups
            .iter()
            .map(|group| column(&format!("ubiquity_{group}")))
            .collect::<Vec<_>>();
        let mut draft = Self {
            member_of: vec![Vec::new(); groups.len()],
            rates: vec![Vec::new(); groups.len()],
            groups,
            ..Self::default()
        };
        for line in lines {
            let fields = line.split('\t').collect::<Vec<_>>();
            let Some(name) = fields.get(name_at) else {
                continue;
            };
            if draft.ids.contains_key(*name) {
                continue;
            }
            let id = draft.names.len() as u16;
            draft.ids.insert((*name).to_string(), id);
            draft.names.push((*name).to_string());
            draft.cutoffs.push(None);
            if let Some(accession) = accession_at.and_then(|at| fields.get(at)) {
                let bare = accession.split('.').next().unwrap_or(accession);
                draft.accessions.insert(bare.to_string(), id);
            }
            for (group, held) in draft.groups.iter().enumerate() {
                let member = set_at
                    .and_then(|at| fields.get(at))
                    .is_some_and(|listed| listed.split(',').any(|entry| entry == held));
                draft.member_of[group].push(member);
                let rate = rate_at[group]
                    .and_then(|at| fields.get(at))
                    .and_then(|field| field.parse::<f64>().ok())
                    .unwrap_or(f64::from(u8::from(member)));
                draft.rates[group].push(rate);
            }
        }
        draft.from_checkm = draft.names.len()..draft.names.len();
        draft
    }

    fn gtdb(&self, model: &str) -> Option<u16> {
        self.ids
            .get(model)
            .copied()
            .filter(|id| (*id as usize) < self.from_checkm.start)
    }

    fn push(&mut self, name: String, cutoff: Option<f32>) -> u16 {
        let id = self.names.len() as u16;
        self.ids.insert(name.clone(), id);
        self.names.push(name);
        self.cutoffs.push(cutoff);
        for (member, rate) in self.member_of.iter_mut().zip(&mut self.rates) {
            member.push(false);
            rate.push(0.0);
        }
        id
    }

    // A CheckM model searched under a GTDB name, or of a family GTDB already holds, is that GTDB
    // marker, so one gene is never two markers of one set.
    fn add_checkm(&mut self, panel: &checkm::Panel, ubiquity: &str) {
        self.families = vec![None; panel.len()];
        for (lineage, model) in panel.members() {
            let Some(group) = self.groups.iter().position(|held| held == lineage) else {
                continue;
            };
            let family = panel.family(model);
            let shared = panel
                .searched_name(model)
                .and_then(|searched| self.gtdb(searched))
                .or_else(|| self.accessions.get(family).copied());
            let marker = match shared {
                Some(marker) => marker,
                None => {
                    let name = format!("{FAMILY_PREFIX}{family}");
                    let marker = match self.ids.get(&name) {
                        Some(marker) => *marker,
                        None => self.push(name, None),
                    };
                    self.families[model as usize] = Some(marker);
                    marker
                }
            };
            self.member_of[group][marker as usize] = true;
        }
        self.from_checkm.end = self.names.len();
        let measured = rates(ubiquity);
        for marker in self.from_checkm.clone() {
            for (group, name) in self.groups.iter().enumerate() {
                let member = self.member_of[group][marker];
                self.rates[group][marker] = measured
                    .get(&(self.names[marker].as_str(), name.as_str()))
                    .copied()
                    .unwrap_or(f64::from(u8::from(member)));
            }
        }
    }

    fn add_panel(&mut self, panel: Option<&Panel>) {
        self.prokaryotic = Some(self.groups.len());
        for lineage in panel.iter().flat_map(|panel| &panel.lineages) {
            let group = self.groups.len();
            self.groups.push(lineage.name.to_string());
            self.member_of.push(vec![false; self.names.len()]);
            self.rates.push(vec![0.0; self.names.len()]);
            for (model, cutoff) in &lineage.models {
                if self.ids.contains_key(model) {
                    continue;
                }
                let id = self.push(model.clone(), Some(*cutoff as f32)) as usize;
                self.member_of[group][id] = true;
                self.rates[group][id] = 1.0;
            }
        }
    }

    fn finish(self) -> MarkerSet {
        let prokaryotic = self.prokaryotic.unwrap_or(self.groups.len());
        let eukaryotic = (0..self.groups.len())
            .map(|group| group >= prokaryotic)
            .collect();
        MarkerSet {
            ids: self.ids,
            names: self.names,
            sets: sets::Sets::new(self.groups, self.member_of, &self.rates)
                .with_eukaryotic(eukaryotic),
            checkm: checkm::Panel::default(),
            cutoffs: self.cutoffs,
            families: self.families,
            from_checkm: self.from_checkm,
        }
    }
}

fn rates(table: &str) -> HashMap<(&str, &str), f64> {
    let mut lines = table.lines();
    let header = lines
        .next()
        .unwrap_or_default()
        .split('\t')
        .collect::<Vec<_>>();
    let mut found = HashMap::new();
    for line in lines {
        let fields = line.split('\t').collect::<Vec<_>>();
        let Some(family) = fields.first() else {
            continue;
        };
        for (column, field) in header.iter().zip(&fields).skip(1) {
            let (Some(group), Ok(rate)) = (column.strip_prefix("ubiquity_"), field.parse()) else {
                continue;
            };
            found.insert((*family, group), rate);
        }
    }
    found
}

impl MarkerSet {
    pub fn embedded() -> Self {
        Self::joined(None)
    }

    pub fn joined(panel: Option<&Panel>) -> Self {
        let mut draft = Draft::parse(TABLE);
        let checkm = checkm::Panel::embedded(|model| draft.gtdb(model));
        draft.add_checkm(&checkm, UBIQUITY);
        draft.add_panel(panel);
        let mut set = draft.finish().with_bounds(SET_TABLE);
        set.checkm = checkm;
        set
    }

    pub fn parse(table: &str) -> Self {
        Self::parse_with(table, None)
    }

    pub fn parse_with(table: &str, panel: Option<&Panel>) -> Self {
        let mut draft = Draft::parse(table);
        draft.add_panel(panel);
        draft.finish()
    }

    pub fn family(&self, model: u16) -> Option<u16> {
        self.families.get(model as usize).copied().flatten()
    }

    pub fn from_checkm(&self, marker: u16) -> bool {
        self.from_checkm.contains(&(marker as usize))
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
