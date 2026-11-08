use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::Result;
use log::debug;

use super::{Chosen, Data, Lineage, fetch, place};
use crate::external::hmmer_engine;
use crate::markers::called::Called;
use crate::markers::checkm::{Counted, reported, score_groups};
use crate::markers::hmm_table;
use crate::markers::{Annotator, Reading};

// CheckM1's lineage_wf defaults for --unique and --multi.
const LEAST_UNIQUE: usize = 10;
const MOST_MULTI: usize = 10;

fn owners(bins: &[&[usize]]) -> HashMap<usize, Vec<usize>> {
    let mut owners = HashMap::<usize, Vec<usize>>::new();
    for (slot, contigs) in bins.iter().enumerate() {
        for contig in *contigs {
            owners.entry(*contig).or_default().push(slot);
        }
    }
    owners
}

// CheckM1 searches each bin alone, so each bin's rows are read as their own table.
fn by_bin(
    table: &str,
    called: &Called,
    owners: &HashMap<usize, Vec<usize>>,
    slots: usize,
) -> Vec<String> {
    let mut held = vec![String::new(); slots];
    for line in table.lines().filter(|line| !line.starts_with('#')) {
        let protein = line
            .split_whitespace()
            .next()
            .and_then(|id| id.parse().ok());
        let Some((contig, _)) = protein.and_then(|protein| called.home(protein)) else {
            continue;
        };
        for slot in owners.get(&contig).into_iter().flatten() {
            held[*slot].push_str(line);
            held[*slot].push('\n');
        }
    }
    held
}

fn searched_only(table: &str, searched: &BTreeSet<String>) -> String {
    table
        .lines()
        .filter(|line| {
            hmm_table::Columns::new(line)
                .at(hmm_table::DOMAIN_MODEL)
                .is_some_and(|model| searched.contains(model))
        })
        .flat_map(|line| [line, "\n"])
        .collect()
}

fn by_model(counted: &[Counted]) -> BTreeMap<u16, Vec<&Counted>> {
    let mut held = BTreeMap::<u16, Vec<&Counted>>::new();
    for copy in counted {
        held.entry(copy.model).or_default().push(copy);
    }
    held
}

fn model_spans(path: &Path) -> Result<HashMap<String, (u64, u64)>> {
    let mut reader = BufReader::new(File::open(path)?);
    let (mut spans, mut line) = (HashMap::new(), String::new());
    let (mut at, mut start, mut name) = (0u64, 0u64, None);
    loop {
        line.clear();
        let read = reader.read_line(&mut line)?;
        if read == 0 {
            break;
        }
        if let Some(rest) = line.strip_prefix("NAME") {
            name = Some(rest.trim().to_string());
        }
        at += read as u64;
        if line.trim_end() == "//" {
            if let Some(held) = name.take() {
                spans.insert(held, (start, at));
            }
            start = at;
        }
    }
    Ok(spans)
}

fn copy_models(
    source: &Path,
    spans: &HashMap<String, (u64, u64)>,
    models: &BTreeSet<&str>,
    target: &Path,
) -> Result<()> {
    let mut source = File::open(source)?;
    let mut sink = BufWriter::new(File::create(target)?);
    for (from, to) in models.iter().filter_map(|model| spans.get(*model)) {
        source.seek(SeekFrom::Start(*from))?;
        std::io::copy(&mut Read::take(&mut source, to - from), &mut sink)?;
    }
    sink.flush()?;
    Ok(())
}

impl Annotator {
    pub fn checkm_lineages(
        &self,
        bins: &[&[usize]],
        names: &[String],
    ) -> Result<Vec<Option<Lineage>>> {
        let mut out = vec![None; bins.len()];
        let Some(data) = self.checkm_tree.as_deref() else {
            return Ok(out);
        };
        let _timer = crate::timing::scope("checkm lineage");
        let folder = tempfile::tempdir()?;
        let contigs = bins
            .iter()
            .flat_map(|bin| bin.iter().copied())
            .collect::<BTreeSet<_>>();
        let called = Called::over(&self.assembly, names, &contigs, folder.path())?;
        let owners = owners(bins);
        let chosen = data.choose(&called, &owners, bins)?;
        let mut by_set = BTreeMap::<&str, Vec<usize>>::new();
        for (slot, held) in chosen.iter().enumerate() {
            if let Some(held) = held {
                by_set.entry(held.uid.as_str()).or_default().push(slot);
            }
        }
        debug!(
            "CheckM lineage sets: {} bins over {} sets",
            bins.len(),
            by_set.len()
        );
        let lineage = data.folder.join(fetch::LINEAGE);
        let spans = model_spans(&lineage)?;
        for (uid, slots) in by_set {
            let picked = slots
                .iter()
                .filter_map(|slot| chosen[*slot].as_ref())
                .collect::<Vec<_>>();
            let widened = picked
                .iter()
                .map(|held| data.relatives.widened(held.groups.iter().flatten()))
                .collect::<Vec<_>>();
            let models = widened
                .iter()
                .flatten()
                .map(String::as_str)
                .collect::<BTreeSet<_>>();
            let hmm = folder.path().join(format!("{uid}.hmm"));
            copy_models(&lineage, &spans, &models, &hmm)?;
            let proteins = folder.path().join(format!("{uid}.faa"));
            let held = slots
                .iter()
                .flat_map(|slot| bins[*slot].iter().copied())
                .collect();
            let panel = data.panel(uid, models.iter().copied());
            let table = match called.copy_contigs(&held, &proteins)? {
                true => hmmer_engine::search(&hmm, &proteins, panel.floor())?,
                false => String::new(),
            };
            let tables = by_bin(&table, &called, &owners, bins.len());
            for ((slot, held), searched) in slots.iter().zip(picked).zip(&widened) {
                let own = searched_only(&tables[*slot], searched);
                let reported = reported(&own, called.proteins_in(bins[*slot]));
                let counted = panel.counted(&reported, |protein| called.home(protein));
                let tally = by_model(&counted);
                let groups = held
                    .groups
                    .iter()
                    .map(|group| {
                        group
                            .iter()
                            .filter_map(|model| panel.id(model))
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>();
                let copies = |model: u16| tally.get(&model).map_or(0, |copies| copies.len() as u32);
                let reading = Reading::of(
                    groups.iter().map(Vec::len).sum(),
                    groups.len(),
                    groups.iter().flatten().map(|model| copies(*model)),
                    score_groups(&groups, copies),
                );
                out[*slot] = Some(Lineage {
                    name: held.name.clone(),
                    reading,
                });
            }
        }
        Ok(out)
    }
}

impl Data {
    fn choose(
        &self,
        called: &Called,
        owners: &HashMap<usize, Vec<usize>>,
        bins: &[&[usize]],
    ) -> Result<Vec<Option<Chosen>>> {
        let slots = bins.len();
        let table = match &called.proteins {
            Some(proteins) => hmmer_engine::search(
                &self.folder.join(fetch::PHYLO),
                proteins,
                self.phylo.floor(),
            )?,
            None => String::new(),
        };
        let tables = by_bin(&table, called, owners, slots);
        let positions = self
            .profiles
            .iter()
            .enumerate()
            .map(|(at, profile)| (profile.name.as_str(), at))
            .collect::<HashMap<_, _>>();
        let mut picked = vec![Vec::new(); self.profiles.len()];
        let mut forced = vec![false; slots];
        let mut file = called.proteins.as_deref().map(File::open).transpose()?;
        for (slot, contigs) in bins.iter().enumerate() {
            let reported = reported(&tables[slot], called.proteins_in(contigs));
            let counted = self
                .phylo
                .counted(&reported, |protein| called.home(protein));
            let held = by_model(&counted);
            let unique = held.values().filter(|copies| copies.len() == 1).count();
            let multi = held.values().filter(|copies| copies.len() > 1).count();
            forced[slot] = unique < LEAST_UNIQUE || multi > MOST_MULTI;
            let Some(file) = file.as_mut() else {
                continue;
            };
            for (model, copies) in held {
                let (Some(at), [copy]) = (positions.get(self.phylo.name(model)), copies.as_slice())
                else {
                    continue;
                };
                picked[*at].push((slot, called.joined(file, copy)?));
            }
        }
        let rows = place::concatenate(&self.profiles, &picked, slots)?;
        Ok(rows
            .iter()
            .zip(forced)
            .map(|(row, forced)| {
                let leaf = row
                    .as_deref()
                    .and_then(|row| self.references.nearest(row))
                    .and_then(|name| self.tree.leaf(name));
                self.tree.choose(leaf, forced)
            })
            .collect())
    }
}
