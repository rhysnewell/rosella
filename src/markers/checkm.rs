use std::collections::HashMap;

use anyhow::Result;

use crate::markers::hmm_table;

pub(crate) const HMM_GZ: &[u8] = include_bytes!("../../data/checkm_markers.hmm.gz");
const SETS: &str = include_str!("../../data/checkm_sets.tsv");
const CLANS: &str = include_str!("../../data/checkm_clans.tsv");

// The rules below are CheckM1's (Parks et al. 2015), so the report reads like it on its own
// marker sets rather than like a different panel scored the same way.
const PSEUDOGENE_SPAN: f64 = 0.3;
const CPR_PREFIX: &str = "cpr_";

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Copies {
    pub set: u8,
    pub model: u16,
    pub copies: u16,
}

struct Lineage {
    name: String,
    gtdb: bool,
    groups: Vec<Vec<u16>>,
}

struct Clan {
    clan: Option<String>,
    nested: Vec<String>,
}

#[derive(Default)]
pub struct Panel {
    ids: HashMap<String, u16>,
    names: Vec<String>,
    pfams: Vec<Option<String>>,
    lineages: Vec<Lineage>,
    clans: HashMap<String, Clan>,
}

impl Panel {
    pub fn embedded(gtdb: impl Fn(&str) -> Option<u16>) -> Self {
        Self::parse(SETS, CLANS, gtdb)
    }

    pub fn parse(sets: &str, clans: &str, gtdb: impl Fn(&str) -> Option<u16>) -> Self {
        let mut panel = Self::default();
        for line in sets.lines().skip(1) {
            let fields = line.split('\t').collect::<Vec<_>>();
            let [set, source, model, group] = fields[..] else {
                continue;
            };
            let Ok(group) = group.parse::<usize>() else {
                continue;
            };
            let from_gtdb = source == "gtdb";
            let id = match from_gtdb {
                true => gtdb(model),
                false => Some(panel.intern(model)),
            };
            let Some(id) = id else {
                continue;
            };
            let at = panel.lineage(set).unwrap_or_else(|| {
                panel.lineages.push(Lineage {
                    name: set.to_string(),
                    gtdb: from_gtdb,
                    groups: Vec::new(),
                });
                panel.lineages.len() - 1
            });
            let groups = &mut panel.lineages[at].groups;
            if groups.len() <= group {
                groups.resize(group + 1, Vec::new());
            }
            groups[group].push(id);
        }
        for line in clans.lines().skip(1) {
            let mut fields = line.split('\t');
            let (Some(pfam), clan, nested) = (fields.next(), fields.next(), fields.next()) else {
                continue;
            };
            panel.clans.insert(
                pfam.to_string(),
                Clan {
                    clan: clan.filter(|clan| !clan.is_empty()).map(str::to_string),
                    nested: nested
                        .unwrap_or_default()
                        .split(',')
                        .filter(|entry| !entry.is_empty())
                        .map(str::to_string)
                        .collect(),
                },
            );
        }
        panel
    }

    fn intern(&mut self, model: &str) -> u16 {
        if let Some(id) = self.ids.get(model) {
            return *id;
        }
        let id = self.names.len() as u16;
        let bare = model.strip_prefix(CPR_PREFIX).unwrap_or(model);
        self.pfams.push(
            bare.starts_with("PF")
                .then(|| bare.split('.').next().unwrap_or(bare).to_string()),
        );
        self.ids.insert(model.to_string(), id);
        self.names.push(model.to_string());
        id
    }

    pub fn lineage(&self, name: &str) -> Option<usize> {
        self.lineages.iter().position(|held| held.name == name)
    }

    pub fn lineage_name(&self, at: usize) -> &str {
        self.lineages
            .get(at)
            .map(|held| held.name.as_str())
            .unwrap_or_default()
    }

    pub fn reads_gtdb(&self, at: usize) -> bool {
        self.lineages.get(at).is_some_and(|held| held.gtdb)
    }

    pub fn id(&self, model: &str) -> Option<u16> {
        self.ids.get(model).copied()
    }

    pub fn name(&self, model: u16) -> &str {
        self.names
            .get(model as usize)
            .map(String::as_str)
            .unwrap_or_default()
    }

    pub fn score(&self, at: usize, copies: impl Fn(u16) -> u32) -> (f64, f64) {
        let Some(lineage) = self.lineages.get(at) else {
            return (0.0, 0.0);
        };
        let groups = lineage.groups.iter().filter(|group| !group.is_empty());
        let (mut present, mut extra, mut count) = (0.0, 0.0, 0usize);
        for group in groups {
            let size = group.len() as f64;
            present += group.iter().filter(|model| copies(**model) > 0).count() as f64 / size;
            extra += group
                .iter()
                .map(|model| f64::from(copies(*model).saturating_sub(1)))
                .sum::<f64>()
                / size;
            count += 1;
        }
        match count {
            0 => (0.0, 0.0),
            _ => (100.0 * present / count as f64, 100.0 * extra / count as f64),
        }
    }

    fn same_clan(&self, ours: u16, theirs: u16) -> bool {
        let (Some(Some(a)), Some(Some(b))) = (
            self.pfams.get(ours as usize),
            self.pfams.get(theirs as usize),
        ) else {
            return false;
        };
        let clan = |pfam: &str| self.clans.get(pfam).and_then(|held| held.clan.as_deref());
        let nested = self
            .clans
            .get(a)
            .is_some_and(|held| held.nested.iter().any(|entry| entry == b));
        clan(a) == clan(b) && !nested
    }

    pub fn tally(
        &self,
        table: &str,
        cutoffs: &Cutoffs,
        contig_of: impl Fn(usize) -> Option<usize>,
        contigs: usize,
    ) -> Vec<Vec<Copies>> {
        let mut best = HashMap::<(u16, usize), Domain>::new();
        for domain in parse(table, self, cutoffs) {
            let held = best.entry((domain.model, domain.protein)).or_insert(domain);
            if held.score < domain.score {
                *held = domain;
            }
        }
        let mut per_contig = vec![HashMap::<(u8, u16), u16>::new(); contigs];
        for (at, lineage) in self.lineages.iter().enumerate() {
            if lineage.gtdb {
                continue;
            }
            let mut members = vec![false; self.names.len()];
            for model in lineage.groups.iter().flatten() {
                members[*model as usize] = true;
            }
            let mut on_protein = HashMap::<usize, Vec<Domain>>::new();
            for domain in best
                .values()
                .filter(|domain| members[domain.model as usize])
            {
                on_protein.entry(domain.protein).or_default().push(*domain);
            }
            let mut proteins = HashMap::<u16, Vec<usize>>::new();
            for (protein, mut held) in on_protein {
                for kept in self.unclashed(&mut held) {
                    proteins.entry(kept.model).or_default().push(protein);
                }
            }
            for (model, mut held) in proteins {
                held.sort_unstable();
                let mut previous: Option<(usize, usize)> = None;
                for protein in held {
                    let Some(contig) = contig_of(protein) else {
                        continue;
                    };
                    if previous.is_some_and(|(seen, home)| seen + 1 == protein && home == contig) {
                        previous = None;
                        continue;
                    }
                    *per_contig[contig].entry((at as u8, model)).or_default() += 1;
                    previous = Some((protein, contig));
                }
            }
        }
        per_contig
            .into_iter()
            .map(|held| {
                let mut copies = held
                    .into_iter()
                    .map(|((set, model), copies)| Copies { set, model, copies })
                    .collect::<Vec<_>>();
                copies.sort_unstable_by_key(|entry| (entry.set, entry.model));
                copies
            })
            .collect()
    }

    fn unclashed(&self, held: &mut [Domain]) -> Vec<Domain> {
        held.sort_by(|a, b| {
            (a.e_value, a.i_evalue)
                .partial_cmp(&(b.e_value, b.i_evalue))
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.model.cmp(&b.model))
        });
        let mut dropped = vec![false; held.len()];
        for ours in 0..held.len() {
            if dropped[ours] || self.pfams[held[ours].model as usize].is_none() {
                continue;
            }
            for theirs in ours + 1..held.len() {
                let (a, b) = (held[ours], held[theirs]);
                if dropped[theirs] || self.pfams[b.model as usize].is_none() {
                    continue;
                }
                let overlap =
                    (a.from <= b.from && a.to > b.from) || (b.from <= a.from && b.to > a.from);
                if overlap && self.same_clan(a.model, b.model) {
                    dropped[theirs] = true;
                }
            }
        }
        held.iter()
            .zip(dropped)
            .filter_map(|(domain, dropped)| (!dropped).then_some(*domain))
            .collect()
    }
}

#[derive(Clone, Copy, Debug)]
struct Domain {
    protein: usize,
    model: u16,
    e_value: f64,
    i_evalue: f64,
    score: f64,
    from: u32,
    to: u32,
}

#[derive(Default)]
pub struct Cutoffs {
    bars: HashMap<String, (f64, f64)>,
}

impl Cutoffs {
    pub fn read(hmm: &std::path::Path) -> Result<Self> {
        Ok(Self::parse(&std::fs::read_to_string(hmm)?))
    }

    pub fn parse(text: &str) -> Self {
        let mut bars = HashMap::new();
        let mut accession = None;
        let mut held = HashMap::<&str, (f64, f64)>::new();
        for line in text.lines() {
            let mut fields = line.split_whitespace();
            match fields.next() {
                Some("ACC") => accession = fields.next().map(str::to_string),
                Some(kind @ ("GA" | "TC" | "NC")) => {
                    let mut values = fields.filter_map(|v| v.trim_end_matches(';').parse().ok());
                    if let (Some(sequence), Some(domain)) = (values.next(), values.next()) {
                        held.insert(kind, (sequence, domain));
                    }
                }
                Some("//") => {
                    if let Some(name) = accession.take() {
                        let tigr = name.contains("TIGR");
                        let order: &[&str] = match tigr {
                            true => &["NC", "GA", "TC"],
                            false => &["GA", "TC", "NC"],
                        };
                        if let Some(bar) = order.iter().find_map(|kind| held.get(kind)) {
                            bars.insert(name, *bar);
                        }
                    }
                    held.clear();
                }
                _ => {}
            }
        }
        Self { bars }
    }

    pub fn floor(&self) -> String {
        let lowest = self
            .bars
            .values()
            .map(|(sequence, domain)| sequence.min(*domain))
            .fold(f64::INFINITY, f64::min);
        match lowest.is_finite() {
            true => format!("{:.2}", (lowest * 100.0).floor() / 100.0),
            false => "0".to_string(),
        }
    }
}

fn parse<'a>(
    table: &'a str,
    panel: &'a Panel,
    cutoffs: &'a Cutoffs,
) -> impl Iterator<Item = Domain> + 'a {
    hmm_table::rows(table).filter_map(|mut fields| {
        let protein = fields.at(hmm_table::DOMAIN_TARGET)?.parse::<usize>().ok()?;
        let accession = fields.at(hmm_table::DOMAIN_ACCESSION)?;
        let length = fields
            .at(hmm_table::DOMAIN_MODEL_LENGTH)?
            .parse::<f64>()
            .ok()?;
        let e_value = fields
            .at(hmm_table::DOMAIN_SEQUENCE_E_VALUE)?
            .parse()
            .ok()?;
        let sequence = fields
            .at(hmm_table::DOMAIN_SEQUENCE_SCORE)?
            .parse::<f64>()
            .ok()?;
        let i_evalue = fields.at(hmm_table::DOMAIN_I_E_VALUE)?.parse().ok()?;
        let score = fields.at(hmm_table::DOMAIN_SCORE)?.parse::<f64>().ok()?;
        let from = fields.at(hmm_table::DOMAIN_ALI_FROM)?.parse::<u32>().ok()?;
        let to = fields.at(hmm_table::DOMAIN_ALI_TO)?.parse::<u32>().ok()?;
        let model = panel.id(accession)?;
        let (bar_sequence, bar_domain) = cutoffs.bars.get(accession)?;
        let aligned = f64::from(to.saturating_sub(from)) / length;
        (aligned >= PSEUDOGENE_SPAN && sequence >= *bar_sequence && score >= *bar_domain).then_some(
            Domain {
                protein,
                model,
                e_value,
                i_evalue,
                score,
                from,
                to,
            },
        )
    })
}
