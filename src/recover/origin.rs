use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::PathBuf;

use anyhow::Result;
use log::{debug, warn};
use needletail::parse_fastx_file;

use crate::cli::{PartitionProbes, RecoverArgs};

pub const BIN_PREFIX: &str = "rosella_bin_";

pub enum Origin<'a> {
    Partition(&'a RecoverArgs),
    Given(Given),
    Ensemble(&'a RecoverArgs, Vec<Given>),
}

pub struct Given {
    names: Vec<String>,
    members: Vec<Vec<String>>,
    prefix: String,
    contamination: HashMap<String, f64>,
    budget: Option<f64>,
}

// Contigs in no given bin are never read or written, so a few bins can be refined alone.
pub struct Placed {
    pub bins: BTreeMap<usize, Vec<usize>>,
    pub contamination: HashMap<usize, f64>,
    pub universe: HashSet<String>,
}

impl Origin<'_> {
    pub fn partition_seeds(&self) -> usize {
        match self {
            Self::Partition(args) | Self::Ensemble(args, _) => args.partition_seeds as usize,
            Self::Given(_) => 1,
        }
    }

    pub fn split_contamination(&self) -> Option<f64> {
        match self {
            Self::Given(given) => given.budget,
            _ => None,
        }
    }

    pub fn probes(&self) -> PartitionProbes {
        match self {
            Self::Partition(args) | Self::Ensemble(args, _) => args.probes.clone(),
            Self::Given(_) => PartitionProbes::default(),
        }
    }

    pub fn prefix(&self) -> String {
        match self {
            Self::Given(given) => given.prefix.clone(),
            _ => BIN_PREFIX.to_string(),
        }
    }

    pub fn universe(&self) -> Option<HashSet<String>> {
        match self {
            Self::Given(given) => Some(given.members.iter().flatten().cloned().collect()),
            _ => None,
        }
    }

    pub fn place(&self, rows: &[String], universe: Option<HashSet<String>>) -> Option<Placed> {
        match (self, universe) {
            (Self::Given(given), Some(universe)) => Some(given.place(rows, universe)),
            _ => None,
        }
    }

    // Each binner's bins are one more arm for the per-bin combination, on the rows the tables hold.
    pub fn candidates(&self, rows: &[String]) -> Vec<BTreeMap<usize, Vec<usize>>> {
        let Self::Ensemble(_, sets) = self else {
            return Vec::new();
        };
        let row_of = row_index(rows);
        sets.iter().map(|set| set.bins_at(&row_of)).collect()
    }
}

fn row_index(rows: &[String]) -> HashMap<&str, usize> {
    rows.iter()
        .enumerate()
        .map(|(row, name)| (name.as_str(), row))
        .collect()
}

impl Given {
    // A contig named in two bins stays in the first, since a run places every contig once.
    pub fn read(
        paths: &[PathBuf],
        prefix: String,
        contamination: HashMap<String, f64>,
        budget: Option<f64>,
    ) -> Result<Self> {
        let mut seen = HashSet::new();
        let mut twice = 0;
        let mut names = Vec::with_capacity(paths.len());
        let mut members = Vec::with_capacity(paths.len());
        for path in paths {
            let mut reader = parse_fastx_file(path)?;
            let mut held = Vec::new();
            while let Some(record) = reader.next() {
                let record = record?;
                let name = crate::contig_id(record.id())?.to_string();
                match seen.insert(name.clone()) {
                    true => held.push(name),
                    false => twice += 1,
                }
            }
            names.push(crate::bins::stem(path));
            members.push(held);
        }
        if twice > 0 {
            warn!("Contigs in two bins: {twice} (first bin kept)");
        }
        Ok(Self {
            names,
            members,
            prefix,
            contamination,
            budget,
        })
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    pub fn place(&self, rows: &[String], universe: HashSet<String>) -> Placed {
        let bins = self.bins_at(&row_index(rows));
        let contamination = self
            .names
            .iter()
            .enumerate()
            .filter_map(|(bin, name)| {
                let found = self.contamination.get(name).copied();
                if found.is_none() && !self.contamination.is_empty() {
                    debug!("Bin quality: no row for {name}");
                }
                Some((bin, found?))
            })
            .collect();
        Placed {
            bins,
            contamination,
            universe,
        }
    }

    fn bins_at(&self, row_of: &HashMap<&str, usize>) -> BTreeMap<usize, Vec<usize>> {
        let mut bins = BTreeMap::new();
        for (bin, held) in self.members.iter().enumerate() {
            let mut contigs = held
                .iter()
                .filter_map(|name| row_of.get(name.as_str()).copied())
                .collect::<Vec<_>>();
            contigs.sort_unstable();
            if !contigs.is_empty() {
                bins.insert(bin, contigs);
            }
        }
        bins
    }
}
