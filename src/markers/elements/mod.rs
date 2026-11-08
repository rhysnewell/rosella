use std::collections::{BTreeSet, HashMap, HashSet};
use std::io::Read;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use log::{debug, warn};

use crate::external::hmmer_engine;
use crate::markers::Annotator;
use crate::markers::called::Called;
use crate::markers::hmm_table::{self, Hits};

mod cache;
mod fetch;

const TABLE_GZ: &[u8] = include_bytes!("../../../data/element_families.tsv.gz");

// The table's own bar on a log likelihood ratio (Antipov et al. 2020).
const BAR: f64 = 7.0;
const CHROMOSOMAL: f64 = std::f64::consts::LN_10;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Call {
    pub virus: f64,
    pub plasmid: f64,
    pub chromosomal: u32,
}

impl Call {
    // A host contig carrying a conjugative element passes the plasmid ratio, so one host family vetoes it.
    pub fn element(&self) -> bool {
        self.virus > BAR || (self.plasmid > BAR && self.chromosomal == 0)
    }
}

struct Weight {
    virus: f64,
    plasmid: f64,
    chromosomal: bool,
}

pub struct Families(HashMap<String, Weight>);

impl Families {
    pub fn parse(table: &str) -> Result<Self> {
        let mut held = HashMap::new();
        for line in table.lines().skip(1).filter(|line| !line.is_empty()) {
            let mut fields = line.split('\t');
            let name = fields.next().context("an empty element family row")?;
            let values = fields
                .map(str::parse::<f64>)
                .collect::<Result<Vec<_>, _>>()
                .with_context(|| format!("element family {name}"))?;
            let [plasmid, chromosome, virus, joined] = values[..] else {
                bail!("element family {name} holds {} rates, not 4", values.len());
            };
            held.insert(
                name.to_string(),
                Weight {
                    virus: virus.ln() - joined.ln(),
                    plasmid: plasmid.ln() - chromosome.ln(),
                    chromosomal: chromosome.ln() - plasmid.max(virus).ln() > CHROMOSOMAL,
                },
            );
        }
        Ok(Self(held))
    }

    fn shipped() -> Result<Self> {
        let mut text = String::new();
        flate2::read::GzDecoder::new(TABLE_GZ).read_to_string(&mut text)?;
        Self::parse(&text)
    }

    fn holds(&self, name: &str) -> bool {
        self.0.contains_key(name)
    }

    // The rates count a family once per sequence, so a repeated family adds nothing.
    pub fn call<'a>(&self, names: impl IntoIterator<Item = &'a str>) -> Call {
        let mut call = Call::default();
        let unique = names.into_iter().collect::<HashSet<_>>();
        for weight in unique.into_iter().filter_map(|name| self.0.get(name)) {
            call.virus += weight.virus;
            call.plasmid += weight.plasmid;
            call.chromosomal += u32::from(weight.chromosomal);
        }
        call
    }
}

pub struct Elements {
    models: PathBuf,
    bars: HashMap<String, (f64, f64)>,
    floor: f64,
    families: Families,
}

impl Elements {
    pub fn prepare(directory: &Path) -> Result<Self> {
        let families = Families::shipped()?;
        let folder = directory.join(fetch::RELEASE);
        if !fetch::complete(&folder) {
            fetch::fetch(directory, &families)?;
        }
        let bars = fetch::read_bars(&folder)?;
        let floor = hmm_table::floor(bars.values().copied(), 1.0);
        Ok(Self {
            models: folder.join(fetch::MODELS),
            bars,
            floor,
            families,
        })
    }

    fn fingerprint() -> u64 {
        crate::digest::fold(TABLE_GZ) ^ crate::digest::fold(fetch::RELEASE.as_bytes())
    }

    fn search(
        &self,
        assembly: &str,
        names: &[String],
        contigs: &BTreeSet<usize>,
    ) -> Result<HashMap<String, Call>> {
        let folder = tempfile::tempdir()?;
        let called = Called::over(assembly, names, contigs, folder.path())?;
        let table = match &called.proteins {
            Some(proteins) => hmmer_engine::search(&self.models, proteins, self.floor)?,
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
                (names[*contig].clone(), self.families.call(held))
            })
            .collect())
    }
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
    pub fn elements(&self, names: &[String], wanted: &BTreeSet<usize>) -> Result<Vec<bool>> {
        let Some(elements) = self.elements.as_deref() else {
            return Ok(Vec::new());
        };
        let _timer = crate::timing::scope("element genes");
        let key = cache::key(self, Elements::fingerprint())?;
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
            calls.extend(elements.search(&self.assembly, names, &missing)?);
            if let Some(directory) = self.cache.as_deref()
                && let Err(error) = cache::write(directory, &key, &calls)
            {
                warn!("Element genes cache: {error}");
            }
        }
        let mut flags = vec![false; names.len()];
        for contig in wanted {
            flags[*contig] = calls.get(&names[*contig]).is_some_and(Call::element);
        }
        debug!(
            "Element genes: {} of {} contigs, {} searched",
            flags.iter().filter(|flag| **flag).count(),
            wanted.len(),
            missing.len()
        );
        Ok(flags)
    }
}
