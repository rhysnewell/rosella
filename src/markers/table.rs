use std::collections::{HashMap, HashSet};

use crate::markers::panels::Panel;
use crate::markers::{checkm, sets};

const TABLE: &str = include_str!("../../data/gtdb_markers.tsv");
const SET_TABLE: &str = include_str!("../../data/marker_sets.tsv");

pub struct MarkerSet {
    names: Names,
    accessions: HashSet<String>,
    pub(crate) sets: sets::Sets,
    pub(crate) checkm: checkm::Panel,
    cutoffs: Vec<Option<f32>>,
}

#[derive(Default)]
pub(super) struct Names {
    ids: HashMap<String, u16>,
    held: Vec<String>,
}

impl Names {
    pub(super) fn push(&mut self, name: &str) -> u16 {
        let id = self.held.len() as u16;
        self.ids.insert(name.to_string(), id);
        self.held.push(name.to_string());
        id
    }

    pub(super) fn id(&self, name: &str) -> Option<u16> {
        self.ids.get(name).copied()
    }

    pub(super) fn name(&self, id: u16) -> &str {
        self.held
            .get(id as usize)
            .map(String::as_str)
            .unwrap_or_default()
    }

    pub(super) fn len(&self) -> usize {
        self.held.len()
    }
}

struct Table<'a> {
    header: Vec<&'a str>,
    rows: Vec<Vec<&'a str>>,
}

impl<'a> Table<'a> {
    fn new(text: &'a str) -> Self {
        let mut lines = text.lines();
        let header = lines.next().unwrap_or_default().split('\t').collect();
        let rows = lines.map(|line| line.split('\t').collect()).collect();
        Self { header, rows }
    }

    fn column(&self, name: &str) -> Option<usize> {
        self.header.iter().position(|field| *field == name)
    }
}

impl MarkerSet {
    pub fn joined(panel: Option<&Panel>) -> Self {
        let mut set = Self::parse_with(TABLE, panel).with_bounds(SET_TABLE);
        set.checkm = checkm::Panel::embedded(|model| set.id(model));
        set
    }

    pub fn parse_with(table: &str, panel: Option<&Panel>) -> Self {
        let table = Table::new(table);
        let Some(name_at) = table.column("model_name") else {
            return Self {
                names: Names::default(),
                accessions: HashSet::new(),
                sets: sets::Sets::default(),
                checkm: checkm::Panel::default(),
                cutoffs: Vec::new(),
            };
        };
        let (set_at, accession_at) = (table.column("sets"), table.column("accession"));
        let mut group_names = Vec::<String>::new();
        let listed = table.rows.iter().filter_map(|fields| fields.get(set_at?));
        for name in listed
            .flat_map(|listed| listed.split(','))
            .filter(|name| !name.is_empty())
        {
            if !group_names.iter().any(|held| held == name) {
                group_names.push(name.to_string());
            }
        }
        let rate_at = group_names
            .iter()
            .map(|group| table.column(&format!("ubiquity_{group}")))
            .collect::<Vec<_>>();

        let mut names = Names::default();
        let mut accessions = HashSet::new();
        let mut member_of = vec![Vec::new(); group_names.len()];
        let mut rates = vec![Vec::new(); group_names.len()];
        for fields in &table.rows {
            let Some(name) = fields.get(name_at) else {
                continue;
            };
            if names.id(name).is_some() {
                continue;
            }
            names.push(name);
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
                if names.id(model).is_some() {
                    continue;
                }
                names.push(model);
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
        let table = Table::new(table);
        let (Some(set_at), Some(max_at)) = (table.column("set"), table.column("max_genome_bp"))
        else {
            return self;
        };
        for fields in &table.rows {
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
        self.names.len() == 0
    }

    pub fn id(&self, model: &str) -> Option<u16> {
        self.names.id(model)
    }

    pub fn name(&self, marker: u16) -> &str {
        self.names.name(marker)
    }
}
