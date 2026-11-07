use std::collections::HashMap;
use std::hash::Hasher;

use crate::markers::hmm_table;
use crate::markers::table::Names;

pub(crate) const HMM_GZ: &[u8] = include_bytes!("../../data/checkm_markers.hmm.gz");
const SETS: &str = include_str!("../../data/checkm_sets.tsv");
const MODELS: &str = include_str!("../../data/checkm_models.tsv");
const CLANS: &str = include_str!("../../data/checkm_clans.tsv");

// The rules below are CheckM1's (Parks et al. 2015), so the report reads like it on its own
// marker sets rather than like a different panel scored the same way.
const PSEUDOGENE_SPAN: f64 = 0.3;
const CPR_PREFIX: &str = "cpr_";

// A CPR tuned model and an older release of one Pfam are the same gene family.
pub fn bare(model: &str) -> &str {
    let model = model.strip_prefix(CPR_PREFIX).unwrap_or(model);
    model.split('.').next().unwrap_or(model)
}

pub(crate) fn fingerprint() -> u64 {
    let mut hash = crate::digest::Fnv1a::default();
    for part in [HMM_GZ, SETS.as_bytes(), MODELS.as_bytes(), CLANS.as_bytes()] {
        hash.write(part);
    }
    hash.finish()
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Copies {
    pub set: u8,
    pub model: u16,
    pub score: f32,
    pub partial: bool,
}

// CheckM1 joins a gene split over two neighbouring calls into one copy, so a copy can hold two.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Counted {
    pub contig: usize,
    pub set: u8,
    pub model: u16,
    pub protein: usize,
    pub partner: Option<usize>,
    pub score: f64,
    pub partial: bool,
}

struct Lineage {
    name: String,
    gtdb: bool,
    groups: Vec<Vec<u16>>,
    models: Vec<u16>,
}

struct Clan {
    clan: Option<String>,
    nested: Vec<String>,
}

#[derive(Default)]
pub struct Panel {
    names: Names,
    bars: Vec<(f64, f64)>,
    searched: HashMap<String, Vec<u16>>,
    searched_as: Vec<Option<String>>,
    pfams: Vec<Option<String>>,
    lineages: Vec<Lineage>,
    clans: HashMap<String, Clan>,
}

impl Panel {
    pub fn embedded(gtdb: impl Fn(&str) -> Option<u16>) -> Self {
        Self::parse(SETS, MODELS, CLANS, gtdb)
    }

    // A CheckM model that is already in the search under another name is not searched twice.
    // `models` maps each one to the name its hits come back under, with its own cutoffs.
    pub fn parse(
        sets: &str,
        models: &str,
        clans: &str,
        gtdb: impl Fn(&str) -> Option<u16>,
    ) -> Self {
        let mut panel = Self::default();
        for line in models.lines().skip(1) {
            let fields = line.split('\t').collect::<Vec<_>>();
            let [model, searched, sequence, domain] = fields[..] else {
                continue;
            };
            let (Ok(sequence), Ok(domain)) = (sequence.parse(), domain.parse()) else {
                continue;
            };
            let id = panel.intern(model);
            panel.bars[id as usize] = (sequence, domain);
            panel.searched_as[id as usize] = Some(searched.to_string());
            panel
                .searched
                .entry(searched.to_string())
                .or_default()
                .push(id);
        }
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
                    models: Vec::new(),
                });
                panel.lineages.len() - 1
            });
            let groups = &mut panel.lineages[at].groups;
            if groups.len() <= group {
                groups.resize(group + 1, Vec::new());
            }
            groups[group].push(id);
        }
        for lineage in &mut panel.lineages {
            lineage.models = lineage.groups.concat();
            lineage.models.sort_unstable();
            lineage.models.dedup();
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
        if let Some(id) = self.names.id(model) {
            return id;
        }
        let family = bare(model);
        self.pfams
            .push(family.starts_with("PF").then(|| family.to_string()));
        self.bars.push((f64::INFINITY, f64::INFINITY));
        self.searched_as.push(None);
        self.names.push(model)
    }

    pub fn floor(&self) -> f64 {
        hmm_table::floor(self.bars.iter().copied(), 1.0)
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

    pub fn models(&self, at: usize) -> &[u16] {
        self.lineages
            .get(at)
            .map_or(&[][..], |held| held.models.as_slice())
    }

    pub fn group_count(&self, at: usize) -> usize {
        self.lineages.get(at).map_or(0, |held| {
            held.groups.iter().filter(|group| !group.is_empty()).count()
        })
    }

    pub fn family(&self, model: u16) -> &str {
        bare(self.name(model))
    }

    pub fn searched_name(&self, model: u16) -> Option<&str> {
        self.searched_as.get(model as usize)?.as_deref()
    }

    pub fn id(&self, model: &str) -> Option<u16> {
        self.names.id(model)
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.len() == 0
    }

    pub fn name(&self, model: u16) -> &str {
        self.names.name(model)
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

    pub fn tally(&self, counted: &[Counted], contigs: usize) -> Vec<Vec<Copies>> {
        let mut per_contig = vec![Vec::new(); contigs];
        for copy in counted {
            per_contig[copy.contig].push(Copies {
                set: copy.set,
                model: copy.model,
                score: copy.score as f32,
                partial: copy.partial,
            });
        }
        for held in &mut per_contig {
            held.sort_unstable_by(in_order);
        }
        per_contig
    }

    pub fn counted(
        &self,
        table: &str,
        contig_and_cut: impl Fn(usize) -> Option<(usize, bool)>,
    ) -> Vec<Counted> {
        let mut best = HashMap::<(u16, usize), Domain>::new();
        for domain in parse(table, self) {
            let held = best.entry((domain.model, domain.protein)).or_insert(domain);
            if held.score < domain.score {
                *held = domain;
            }
        }
        let mut counted: Vec<Counted> = Vec::new();
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
            let mut proteins = HashMap::<u16, Vec<(usize, f64)>>::new();
            for (protein, mut held) in on_protein {
                for kept in self.unclashed(&mut held) {
                    proteins
                        .entry(kept.model)
                        .or_default()
                        .push((protein, kept.score));
                }
            }
            for (model, mut held) in proteins {
                held.sort_unstable_by_key(|(protein, _)| *protein);
                let mut previous: Option<(usize, usize)> = None;
                for (protein, score) in held {
                    let Some((contig, cut)) = contig_and_cut(protein) else {
                        continue;
                    };
                    if previous.is_some_and(|(seen, home)| seen + 1 == protein && home == contig) {
                        if let Some(last) = counted.last_mut() {
                            last.partner = Some(protein);
                            last.score = last.score.max(score);
                            last.partial |= cut;
                        }
                        previous = None;
                        continue;
                    }
                    counted.push(Counted {
                        contig,
                        set: at as u8,
                        model,
                        protein,
                        partner: None,
                        score,
                        partial: cut,
                    });
                    previous = Some((protein, contig));
                }
            }
        }
        counted
    }

    fn unclashed(&self, held: &mut [Domain]) -> Vec<Domain> {
        held.sort_by(|a, b| {
            a.e_value
                .total_cmp(&b.e_value)
                .then(a.i_evalue.total_cmp(&b.i_evalue))
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

// A contig's copies of one model sit together with the best first, which is the copy the bin reads.
fn in_order(a: &Copies, b: &Copies) -> std::cmp::Ordering {
    (a.set, a.model)
        .cmp(&(b.set, b.model))
        .then(b.score.total_cmp(&a.score))
        .then(a.partial.cmp(&b.partial))
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

fn parse<'a>(table: &'a str, panel: &'a Panel) -> impl Iterator<Item = Domain> + 'a {
    hmm_table::domains(table).flat_map(|row| {
        let (Some(ids), Some(e_value), Some(i_evalue)) = (
            panel.searched.get(row.model),
            row.sequence_e_value,
            row.i_e_value,
        ) else {
            return Vec::new();
        };
        let (from, to) = (row.reach.protein_from, row.reach.protein_to);
        let aligned = f64::from(to.saturating_sub(from)) / f64::from(row.reach.model_length);
        ids.iter()
            .filter(|id| {
                let (bar_sequence, bar_domain) = panel.bars[**id as usize];
                aligned >= PSEUDOGENE_SPAN
                    && row.sequence_score >= bar_sequence
                    && row.score >= bar_domain
            })
            .map(|id| Domain {
                protein: row.protein,
                model: *id,
                e_value,
                i_evalue,
                score: row.score,
                from,
                to,
            })
            .collect::<Vec<_>>()
    })
}
