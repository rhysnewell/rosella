use std::collections::{BTreeSet, HashMap, HashSet};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use log::debug;

use super::nodes::parse_newick;
use crate::external::hmmer_engine;
use crate::markers::fragments::cutoffs;

pub const RELEASE: &str = "checkm_data_2015_01_16";
// Every node's models are kept as clan relatives, so a cache pruned to the selected sets is refetched.
pub const FOLDER: &str = "checkm_data_2015_01_16_paths";
const SOURCE: &str = "https://data.ace.uq.edu.au/public/CheckM_databases";
pub const SELECTED: &str = include_str!("../../../data/checkm_selected.tsv");

pub const PHYLO: &str = "phylo.hmm";
pub const LINEAGE: &str = "lineage.hmm";
pub const TREE: &str = "genome_tree.tre";
pub const REFERENCES: &str = "genome_tree.fasta";
pub const NODES: &str = "nodes.tsv";
pub const MODELS: &str = "models.tsv";
pub const CLANS: &str = "clans.tsv";
const FILES: [&str; 7] = [PHYLO, LINEAGE, TREE, REFERENCES, NODES, MODELS, CLANS];

const METADATA_ENTRY: &str = "genome_tree/genome_tree.metadata.tsv";
const MISSING_ENTRY: &str = "genome_tree/missing_duplicate_genes_50.tsv";
const TREE_ENTRY: &str = "genome_tree/genome_tree_full.refpkg/genome_tree.tre";
const REFERENCES_ENTRY: &str = "genome_tree/genome_tree_full.refpkg/genome_tree.fasta";
const PHYLO_ENTRY: &str = "hmms/phylo.hmm";
const MODELS_ENTRY: &str = "hmms/checkm.hmm";
const PFAM_ENTRY: &str = "pfam/Pfam-A.hmm.dat";

pub fn complete(folder: &Path) -> bool {
    FILES.iter().all(|file| folder.join(file).exists())
}

pub fn fetch(directory: &Path) -> Result<()> {
    std::fs::create_dir_all(directory)?;
    let url = format!("{SOURCE}/{RELEASE}.tar.gz");
    debug!("CheckM lineage sets: fetching {url}");
    let staging = tempfile::tempdir_in(directory)?;
    let folder = staging.path().join(RELEASE);
    std::fs::create_dir(&folder)?;
    let body = ureq::get(&url)
        .call()
        .with_context(|| format!("could not fetch {url}"))?
        .into_body();
    let mut archive = tar::Archive::new(flate2::read::GzDecoder::new(body.into_reader()));
    let mut pruning = Pruning::default();
    for entry in archive.entries()? {
        let mut entry = entry?;
        let path = entry.path()?.to_string_lossy().into_owned();
        pruning.take(path.trim_start_matches("./"), &mut entry, &folder)?;
    }
    pruning.finish(&folder)?;
    let target = directory.join(FOLDER);
    match std::fs::rename(&folder, &target) {
        Ok(()) => Ok(()),
        Err(_) if complete(&target) => Ok(()),
        Err(error) => Err(error).with_context(|| format!("could not fill {}", target.display())),
    }
}

struct Row {
    genomes: u32,
    taxonomy: String,
    sets: String,
}

#[derive(Default)]
struct Pruning {
    metadata: HashMap<String, Row>,
    removed: HashMap<String, Vec<String>>,
    uids: Vec<String>,
    candidates: HashMap<String, Vec<Vec<String>>>,
    models: Vec<(String, (f64, f64))>,
    clans: String,
    seen: HashSet<&'static str>,
}

fn quoted(text: &str) -> impl Iterator<Item = &str> {
    text.split('\'').skip(1).step_by(2)
}

fn groups(sets: &str) -> Vec<Vec<String>> {
    sets.split("set(")
        .skip(1)
        .map(|group| quoted(group).map(str::to_string).collect())
        .collect()
}

fn header<'a>(block: &'a str, key: &str) -> Option<&'a str> {
    block.lines().find_map(|line| {
        let rest = line.strip_prefix(key)?;
        rest.starts_with(char::is_whitespace).then(|| rest.trim())
    })
}

fn bars(block: &str, accession: &str) -> Option<(f64, f64)> {
    let noise = header(block, "NC").and_then(cutoffs);
    let gathering = header(block, "GA").and_then(cutoffs);
    let trusted = header(block, "TC").and_then(cutoffs);
    match accession.contains("TIGR") {
        true => noise.or(gathering).or(trusted),
        false => gathering.or(trusted).or(noise),
    }
}

// Named by accession, so a table row names the model the way the sets do.
fn renamed(block: &str) -> Option<(String, String)> {
    let accession = header(block, "ACC")?.to_string();
    let text = block
        .lines()
        .map(|line| match line.starts_with("NAME") {
            true => format!("NAME  {accession}\n"),
            false => format!("{line}\n"),
        })
        .collect();
    Some((accession, text))
}

impl Pruning {
    fn take(&mut self, path: &str, entry: &mut impl Read, folder: &Path) -> Result<()> {
        let entry_name = match path {
            METADATA_ENTRY => METADATA_ENTRY,
            MISSING_ENTRY => MISSING_ENTRY,
            TREE_ENTRY => TREE_ENTRY,
            REFERENCES_ENTRY => REFERENCES_ENTRY,
            PHYLO_ENTRY => PHYLO_ENTRY,
            MODELS_ENTRY => MODELS_ENTRY,
            PFAM_ENTRY => PFAM_ENTRY,
            _ => return Ok(()),
        };
        self.seen.insert(entry_name);
        match entry_name {
            METADATA_ENTRY => self.read_metadata(&text(entry)?),
            MISSING_ENTRY => {
                for line in text(entry)?.lines() {
                    let mut fields = line.split('\t');
                    let uid = fields.next().unwrap_or_default().to_string();
                    let genes = fields.flat_map(quoted).map(as_model).collect();
                    self.removed.insert(uid, genes);
                }
            }
            TREE_ENTRY => {
                let newick = text(entry)?;
                std::fs::write(folder.join(TREE), &newick)?;
                let topology = parse_newick(&newick)?;
                self.uids = topology
                    .names
                    .iter()
                    .filter_map(|name| name.split_once('|').map(|(uid, _)| uid.to_string()))
                    .collect();
            }
            REFERENCES_ENTRY => {
                let mut sink = BufWriter::new(std::fs::File::create(folder.join(REFERENCES))?);
                std::io::copy(entry, &mut sink)?;
                sink.flush()?;
            }
            PHYLO_ENTRY => self.keep_models(entry, &folder.join(PHYLO), |_| true)?,
            MODELS_ENTRY => {
                self.choose_candidates()?;
                let wanted = self
                    .candidates
                    .values()
                    .flatten()
                    .flatten()
                    .cloned()
                    .collect::<HashSet<_>>();
                self.keep_models(entry, &folder.join(LINEAGE), |model| wanted.contains(model))?;
            }
            _ => self.clans = text(entry)?,
        }
        Ok(())
    }

    fn read_metadata(&mut self, text: &str) {
        for line in text.lines().skip(1) {
            let fields = line.split('\t').collect::<Vec<_>>();
            let (Some(uid), Some(genomes), Some(taxonomy), Some(sets)) =
                (fields.first(), fields.get(1), fields.get(2), fields.get(10))
            else {
                continue;
            };
            let row = Row {
                genomes: genomes.parse().unwrap_or_default(),
                taxonomy: taxonomy.to_string(),
                sets: sets.to_string(),
            };
            self.metadata.insert(uid.to_string(), row);
        }
    }

    fn choose_candidates(&mut self) -> Result<()> {
        if self.uids.is_empty() || self.metadata.is_empty() {
            bail!("the CheckM archive lists its models before its genome tree");
        }
        for uid in &self.uids {
            if let Some(row) = self.metadata.get(uid) {
                self.candidates.insert(uid.to_string(), groups(&row.sets));
            }
        }
        Ok(())
    }

    fn keep_models(
        &mut self,
        entry: &mut impl Read,
        path: &Path,
        wanted: impl Fn(&str) -> bool,
    ) -> Result<()> {
        let mut sink = BufWriter::new(std::fs::File::create(path)?);
        hmmer_engine::each_model(BufReader::new(entry), |_, block| {
            let Some((accession, text)) = renamed(&block) else {
                return Ok(());
            };
            if !wanted(&accession) {
                return Ok(());
            }
            let found = bars(&block, &accession)
                .with_context(|| format!("{accession} carries no cutoff"))?;
            sink.write_all(text.as_bytes())?;
            self.models.push((accession, found));
            Ok(())
        })?;
        sink.flush()?;
        Ok(())
    }

    fn finish(self, folder: &Path) -> Result<()> {
        for wanted in [
            METADATA_ENTRY,
            MISSING_ENTRY,
            TREE_ENTRY,
            REFERENCES_ENTRY,
            PHYLO_ENTRY,
            MODELS_ENTRY,
            PFAM_ENTRY,
        ] {
            if !self.seen.contains(wanted) {
                bail!("the CheckM archive holds no {wanted}");
            }
        }
        let mut nodes = String::from("uid\tgenomes\ttaxonomy\tremoved\tgroups\n");
        for uid in &self.uids {
            let Some(row) = self.metadata.get(uid) else {
                continue;
            };
            let removed = self.removed.get(uid).map(|genes| genes.join(","));
            let groups = self.candidates.get(uid).map(|groups| {
                groups
                    .iter()
                    .map(|group| group.join(","))
                    .collect::<Vec<_>>()
                    .join(";")
            });
            nodes.push_str(&format!(
                "{uid}\t{}\t{}\t{}\t{}\n",
                row.genomes,
                row.taxonomy,
                removed.unwrap_or_default(),
                groups.unwrap_or_default(),
            ));
        }
        std::fs::write(folder.join(NODES), nodes)?;
        let mut models = String::from("model\tsearched\tsequence\tdomain\n");
        for (model, (sequence, domain)) in &self.models {
            models.push_str(&format!("{model}\t{model}\t{sequence}\t{domain}\n"));
        }
        std::fs::write(folder.join(MODELS), models)?;
        let pfams = self
            .models
            .iter()
            .map(|(model, _)| crate::markers::checkm::bare(model))
            .filter(|model| model.starts_with("PF"))
            .collect::<HashSet<_>>();
        std::fs::write(folder.join(CLANS), clans(&self.clans, &pfams))?;
        Ok(())
    }
}

fn text(entry: &mut impl Read) -> Result<String> {
    let mut held = String::new();
    entry.read_to_string(&mut held)?;
    Ok(held)
}

fn as_model(gene: &str) -> String {
    match gene.strip_prefix("pfam") {
        Some(number) => format!("PF{number}"),
        None => gene.to_string(),
    }
}

// Read as CheckM1 reads Pfam-A.hmm.dat: a clan attaches to the last accession seen, and nesting
// runs both ways between the two families.
fn clans(text: &str, wanted: &HashSet<&str>) -> String {
    let (mut accession, mut name) = (String::new(), String::new());
    let mut ids = HashMap::<String, String>::new();
    let mut clan = HashMap::<String, String>::new();
    let mut nested = HashMap::<String, BTreeSet<String>>::new();
    let mut pairs = Vec::new();
    for line in text.lines() {
        let mut fields = line.split_whitespace();
        let (Some("#=GF"), Some(key), Some(value)) = (fields.next(), fields.next(), fields.next())
        else {
            continue;
        };
        match key {
            "ID" => name = value.to_string(),
            "AC" => {
                accession = value.split('.').next().unwrap_or(value).to_string();
                ids.insert(name.clone(), accession.clone());
            }
            "CL" => {
                clan.insert(accession.clone(), value.to_string());
            }
            "NE" => pairs.push((name.clone(), value.to_string())),
            _ => {}
        }
    }
    for (ours, theirs) in pairs {
        let (Some(ours), Some(theirs)) = (ids.get(&ours), ids.get(&theirs)) else {
            continue;
        };
        nested
            .entry(ours.clone())
            .or_default()
            .insert(theirs.clone());
        nested
            .entry(theirs.clone())
            .or_default()
            .insert(ours.clone());
    }
    let mut rows = String::from("pfam\tclan\tnested\n");
    let mut kept = wanted.iter().collect::<Vec<_>>();
    kept.sort_unstable();
    for pfam in kept {
        let held = nested.get(*pfam).map(|held| {
            held.iter()
                .map(String::as_str)
                .collect::<Vec<_>>()
                .join(",")
        });
        rows.push_str(&format!(
            "{pfam}\t{}\t{}\n",
            clan.get(*pfam).map_or("", String::as_str),
            held.unwrap_or_default(),
        ));
    }
    rows
}
