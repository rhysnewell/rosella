use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::BufRead;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use log::{debug, warn};
use needletail::parse_fastx_file;

use crate::cli::RecoverArgs;

pub const BIN_PREFIX: &str = "rosella_bin_";

pub enum Origin<'a> {
    Partition(&'a RecoverArgs),
    Given(Given),
    Ensemble(&'a RecoverArgs, Offered),
}

// With `finished`, the first set is the run's own finished bins and no partition is made.
pub struct Offered {
    pub sets: Vec<(String, Given)>,
    pub finished: bool,
}

pub struct Given {
    names: Vec<String>,
    members: Vec<Vec<String>>,
    prefix: String,
    contamination: HashMap<String, f64>,
    budget: Option<f64>,
}

pub struct BinSet {
    pub name: String,
    pub bins: Vec<(String, Vec<usize>)>,
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

    // Given bins teach a space and call genes on those bins alone, which serves no later run.
    pub fn keeps_cache(&self) -> bool {
        !matches!(self, Self::Given(_))
    }

    pub fn finished(&self) -> bool {
        matches!(self, Self::Ensemble(_, offered) if offered.finished)
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
    pub fn candidates(&self, rows: &[String]) -> Vec<BinSet> {
        let Self::Ensemble(_, offered) = self else {
            return Vec::new();
        };
        let row_of = row_index(rows);
        offered
            .sets
            .iter()
            .map(|(name, set)| BinSet {
                name: name.clone(),
                bins: set
                    .bins_at(&row_of, name)
                    .into_iter()
                    .map(|(bin, contigs)| (set.names[bin].clone(), contigs))
                    .collect(),
            })
            .collect()
    }
}

fn row_index(rows: &[String]) -> HashMap<&str, usize> {
    crate::rows::positions(rows.iter().map(String::as_str))
}

impl Given {
    // A contig named in two bins stays in the first, since a run places every contig once.
    pub fn read(
        paths: &[PathBuf],
        prefix: String,
        contamination: HashMap<String, f64>,
        budget: Option<f64>,
    ) -> Result<Self> {
        let mut bins = Vec::with_capacity(paths.len());
        for path in paths {
            let mut reader = parse_fastx_file(path)?;
            let mut held = Vec::new();
            while let Some(record) = reader.next() {
                held.push(crate::contig_id(record?.id())?.to_string());
            }
            bins.push((crate::bins::stem(path), held));
        }
        Ok(Self::of(bins, prefix, contamination, budget))
    }

    // Lines starting `@` or `#` are skipped, so a CAMI binning reads as a plain contig and bin table.
    pub fn table(path: &Path) -> Result<Self> {
        let mut at = HashMap::new();
        let mut bins: Vec<(String, Vec<String>)> = Vec::new();
        for line in crate::get_file_reader(path)?.lines() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with(['@', '#']) {
                continue;
            }
            let mut fields = line.split('\t');
            let (Some(contig), Some(bin)) = (fields.next(), fields.next()) else {
                bail!(
                    "{}: expected a contig and a bin, tab separated, in {line:?}",
                    path.display()
                );
            };
            let slot = *at.entry(bin.trim().to_string()).or_insert_with(|| {
                bins.push((bin.trim().to_string(), Vec::new()));
                bins.len() - 1
            });
            bins[slot]
                .1
                .push(crate::contig_id(contig.as_bytes())?.to_string());
        }
        Ok(Self::of(bins, String::new(), HashMap::new(), None))
    }

    pub fn offered(source: &str) -> Result<Self> {
        match Path::new(source).is_dir() {
            true => Self::read(
                &crate::bins::fasta_in(source)?,
                String::new(),
                HashMap::new(),
                None,
            ),
            false => Self::table(Path::new(source)),
        }
    }

    // A recover output writes loose contigs to files beside its bins, so its bins are the ones its
    // quality table names.
    pub fn recovered(source: &str) -> Result<Self> {
        let folder = Path::new(source);
        if !folder.is_dir() {
            return Self::table(folder);
        }
        let table = folder.join(crate::defaults::QUALITY_FILE);
        let rows = std::fs::read_to_string(&table).with_context(|| {
            format!(
                "{} is not a rosella recover output: it holds no {}",
                folder.display(),
                crate::defaults::QUALITY_FILE
            )
        })?;
        let files = rows
            .lines()
            .skip(1)
            .filter_map(|row| row.split('\t').next())
            .filter(|bin| !bin.is_empty())
            .map(|bin| folder.join(format!("{bin}.{}", crate::defaults::FASTA_EXTENSION)))
            .collect::<Vec<_>>();
        Self::read(&files, String::new(), HashMap::new(), None)
    }

    fn of(
        bins: Vec<(String, Vec<String>)>,
        prefix: String,
        contamination: HashMap<String, f64>,
        budget: Option<f64>,
    ) -> Self {
        let mut seen = HashSet::new();
        let mut twice = 0;
        let mut names = Vec::with_capacity(bins.len());
        let mut members = Vec::with_capacity(bins.len());
        for (name, held) in bins {
            let before = held.len();
            let held = held
                .into_iter()
                .filter(|contig| seen.insert(contig.clone()))
                .collect::<Vec<_>>();
            twice += before - held.len();
            names.push(name);
            members.push(held);
        }
        if twice > 0 {
            warn!("Contigs in two bins: {twice} (first bin kept)");
        }
        Self {
            names,
            members,
            prefix,
            contamination,
            budget,
        }
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    pub fn place(&self, rows: &[String], universe: HashSet<String>) -> Placed {
        let bins = self.bins_at(&row_index(rows), "Given bins");
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

    fn bins_at(&self, row_of: &HashMap<&str, usize>, source: &str) -> BTreeMap<usize, Vec<usize>> {
        let mut bins = BTreeMap::new();
        let mut missing = 0;
        for (bin, held) in self.members.iter().enumerate() {
            let mut contigs = held
                .iter()
                .filter_map(|name| row_of.get(name.as_str()).copied())
                .collect::<Vec<_>>();
            missing += held.len() - contigs.len();
            contigs.sort_unstable();
            if !contigs.is_empty() {
                bins.insert(bin, contigs);
            }
        }
        if missing > 0 {
            let total = self.members.iter().map(Vec::len).sum::<usize>();
            let lost = format!(
                "{source}: {missing} of {total} binned contigs are not in the run (renamed, blank, \
                 or under the attach floor)"
            );
            // Other binners keep contigs under the floor, so only bins that mostly miss the
            // assembly are worth a warning.
            match 2 * missing > total {
                true => warn!("{lost}"),
                false => debug!("{lost}"),
            }
        }
        bins
    }
}
