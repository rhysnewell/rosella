use std::collections::{BTreeSet, HashMap, HashSet};
use std::io::Read;

use anyhow::{Context, Result, bail};
use log::{debug, warn};

use crate::external::hmmer_engine;
use crate::markers::Annotator;
use crate::markers::called::Called;
use crate::markers::fragments::cutoffs;
use crate::markers::hmm_table::{self, Hits};

mod cache;

// Only phage and plasmid families ship, since CheckM markers and coding density already mark host.
const MODELS_ZST: &[u8] = include_bytes!("../../../data/virus_families.hmm.zst");
const TABLE_GZ: &[u8] = include_bytes!("../../../data/virus_families.tsv.gz");

// The table's own bar on a log likelihood ratio (Antipov et al. 2020).
const BAR: f64 = 7.0;

// Only viruses leave. Gold standards disagree on whether a plasmid belongs to its host genome.
pub fn viral(score: f64) -> bool {
    score > BAR
}

pub struct Families(HashMap<String, f64>);

impl Families {
    pub fn parse(table: &str) -> Result<Self> {
        let mut held = HashMap::new();
        for line in table.lines().skip(1).filter(|line| !line.is_empty()) {
            let mut fields = line.split('\t');
            let name = fields.next().context("an empty virus family row")?;
            let values = fields
                .map(str::parse::<f64>)
                .collect::<Result<Vec<_>, _>>()
                .with_context(|| format!("virus family {name}"))?;
            let [_, _, virus, joined] = values[..] else {
                bail!("virus family {name} holds {} rates, not 4", values.len());
            };
            held.insert(name.to_string(), virus.ln() - joined.ln());
        }
        Ok(Self(held))
    }

    // The rates count a family once per sequence, so a repeated family adds nothing.
    pub fn virus<'a>(&self, names: impl IntoIterator<Item = &'a str>) -> f64 {
        let unique = names.into_iter().collect::<HashSet<_>>();
        unique.into_iter().filter_map(|name| self.0.get(name)).sum()
    }
}

pub struct Viruses {
    models: Vec<u8>,
    bars: HashMap<String, (f64, f64)>,
    floor: f64,
    families: Families,
}

impl Viruses {
    pub fn shipped() -> Result<Self> {
        let models = zstd::decode_all(MODELS_ZST)?;
        let mut table = String::new();
        flate2::read::GzDecoder::new(TABLE_GZ).read_to_string(&mut table)?;
        let bars = gathering(std::str::from_utf8(&models)?);
        let floor = hmm_table::floor(bars.values().copied(), 1.0);
        Ok(Self {
            models,
            bars,
            floor,
            families: Families::parse(&table)?,
        })
    }

    fn fingerprint() -> u64 {
        crate::digest::fold(MODELS_ZST) ^ crate::digest::fold(TABLE_GZ)
    }

    fn search(
        &self,
        assembly: &str,
        names: &[String],
        contigs: &BTreeSet<usize>,
    ) -> Result<HashMap<String, f64>> {
        let folder = tempfile::tempdir()?;
        let called = Called::over(assembly, names, contigs, folder.path())?;
        let hmm = folder.path().join("viruses.hmm");
        std::fs::write(&hmm, &self.models)?;
        let table = match &called.proteins {
            Some(proteins) => hmmer_engine::search(&hmm, proteins, self.floor, "virus families")?,
            None => String::new(),
        };
        let mut found = HashMap::<usize, Vec<String>>::new();
        for (protein, best) in gathered(&table, &self.bars) {
            if let Some((contig, _)) = called.home(protein) {
                found.entry(contig).or_default().push(best.model);
            }
        }
        Ok(contigs
            .iter()
            .map(|contig| {
                let held = found.get(contig).into_iter().flatten().map(String::as_str);
                (names[*contig].clone(), self.families.virus(held))
            })
            .collect())
    }
}

fn gathering(models: &str) -> HashMap<String, (f64, f64)> {
    let mut bars = HashMap::new();
    let mut name = None;
    for line in models.lines() {
        if let Some(rest) = line.strip_prefix("NAME") {
            name = Some(rest.trim().to_string());
        } else if let Some(found) = line.strip_prefix("GA ").and_then(cutoffs)
            && let Some(name) = name.take()
        {
            bars.insert(name, found);
        }
    }
    bars
}

fn gathered(table: &str, bars: &HashMap<String, (f64, f64)>) -> Hits {
    let mut best = Hits::new();
    for domain in hmm_table::domains(table) {
        let Some((sequence, single)) = bars.get(domain.model) else {
            continue;
        };
        if domain.sequence_score >= *sequence && domain.score >= *single {
            hmm_table::keep_best(
                &mut best,
                domain.protein,
                domain.model,
                domain.sequence_score,
                domain.reach,
            );
        }
    }
    best
}

impl Annotator {
    pub fn viral(&self, names: &[String], wanted: &BTreeSet<usize>) -> Result<Vec<bool>> {
        let Some(viruses) = self.viruses.as_deref() else {
            return Ok(Vec::new());
        };
        let _timer = crate::timing::scope("viruses");
        let key = cache::key(self, Viruses::fingerprint())?;
        let mut calls = self
            .cache
            .as_deref()
            .map(|directory| cache::read(directory, &key))
            .unwrap_or_default();
        let missing = wanted
            .iter()
            .copied()
            .filter(|contig| !calls.contains_key(&names[*contig]))
            .collect::<BTreeSet<_>>();
        if !missing.is_empty() {
            calls.extend(viruses.search(&self.assembly, names, &missing)?);
            if let Some(directory) = self.cache.as_deref()
                && let Err(error) = cache::write(directory, &key, &calls)
            {
                warn!("Virus cache: {error}");
            }
        }
        let mut flags = vec![false; names.len()];
        for contig in wanted {
            flags[*contig] = calls
                .get(&names[*contig])
                .is_some_and(|score| viral(*score));
        }
        debug!(
            "Viral contigs: {} of {}, {} searched",
            flags.iter().filter(|flag| **flag).count(),
            wanted.len(),
            missing.len()
        );
        Ok(flags)
    }
}
