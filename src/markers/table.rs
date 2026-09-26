use std::collections::HashMap;

use crate::markers::sets;

const TABLE: &str = include_str!("../../data/gtdb_markers.tsv");
const SET_TABLE: &str = include_str!("../../data/marker_sets.tsv");

pub struct MarkerSet {
    ids: HashMap<String, u16>,
    names: Vec<String>,
    pub(crate) sets: sets::Sets,
}

const FALLBACK_SETS: [(&str, &str); 2] = [("bac", "bac120"), ("ar", "ar53")];

impl MarkerSet {
    pub fn embedded() -> Self {
        Self::parse(TABLE).with_bounds(SET_TABLE)
    }

    pub fn parse(table: &str) -> Self {
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
                sets: sets::Sets::default(),
            };
        };
        let set_at = column("sets");
        let domain_at = column("domain");
        let group_names = match set_at {
            Some(_) => set_names(table),
            None => FALLBACK_SETS
                .iter()
                .map(|(set, _)| (*set).to_string())
                .collect(),
        };
        let rate_at = group_names
            .iter()
            .map(|group| column(&format!("ubiquity_{group}")))
            .collect::<Vec<_>>();
        let copies_at = group_names
            .iter()
            .map(|group| column(&format!("single_copy_{group}")))
            .collect::<Vec<_>>();

        let mut ids = HashMap::new();
        let mut names = Vec::new();
        let mut member_of = vec![Vec::new(); group_names.len()];
        let mut rates = vec![Vec::new(); group_names.len()];
        let mut copies = vec![Vec::new(); group_names.len()];
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
            for (group, held) in group_names.iter().enumerate() {
                let member = match set_at.and_then(|at| fields.get(at)) {
                    Some(listed) => listed.split(',').any(|entry| entry == held),
                    None => domain_at
                        .and_then(|at| fields.get(at))
                        .is_some_and(|domain| domain.contains(FALLBACK_SETS[group].1)),
                };
                member_of[group].push(member);
                let rate = rate_at[group]
                    .and_then(|at| fields.get(at))
                    .and_then(|field| field.parse::<f64>().ok())
                    .unwrap_or(f64::from(u8::from(member)));
                rates[group].push(rate);
                // A table with no single copy column weights every duplicate in full, which
                // is what the scorer did before the column existed.
                copies[group].push(
                    copies_at[group]
                        .and_then(|at| fields.get(at))
                        .and_then(|field| field.parse::<f64>().ok())
                        .unwrap_or(1.0),
                );
            }
        }
        Self {
            ids,
            names,
            sets: sets::Sets::new(group_names, member_of, &rates, copies),
        }
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
