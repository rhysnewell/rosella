use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{Context, Result, bail};
use log::info;

use super::fragments::{self, Bar, Bars};
use super::{Hit, MarkerSet};
use crate::digest::fold;
use crate::external::hmmer_engine;
use crate::quality::orfs::Orf;

// BUSCO's sets are CC BY-ND, so they are fetched as released and never shipped.
const SOURCE: &str = "https://busco-data.ezlab.org/v5/data/lineages";
const RELEASE: &str = "2026-05-22";
const ORTHODB: &str = "odb12";
pub const LINEAGES: [&str; 4] = ["eukaryota", "fungi", "viridiplantae", "metazoa"];
const JOINED: &str = "eukaryotes.hmm";
const CUTOFFS: &str = "scores_cutoff";
const HMMS: &str = "hmms";

pub struct Lineage {
    pub name: &'static str,
    pub models: Vec<(String, f64)>,
}

pub struct Panel {
    hmm: PathBuf,
    pub lineages: Vec<Lineage>,
    fingerprint: u64,
}

impl Panel {
    pub fn prepare(directory: &Path) -> Result<Self> {
        std::fs::create_dir_all(directory)?;
        let hmm = directory.join(JOINED);
        let mut lineages = Vec::new();
        let mut cutoffs = Vec::new();
        for name in LINEAGES {
            let folder = directory.join(format!("{name}_{ORTHODB}"));
            if !folder.join(CUTOFFS).exists() || !(hmm.exists() || folder.join(HMMS).exists()) {
                fetch(directory, name)?;
            }
            let text = std::fs::read_to_string(folder.join(CUTOFFS))?;
            let models = text
                .lines()
                .filter_map(|line| {
                    let mut fields = line.split_whitespace();
                    Some((fields.next()?.to_string(), fields.next()?.parse().ok()?))
                })
                .collect::<Vec<_>>();
            cutoffs.extend_from_slice(text.as_bytes());
            lineages.push(Lineage { name, models });
        }
        if !hmm.exists() {
            join(directory, &lineages, &hmm)?;
            for name in LINEAGES {
                std::fs::remove_dir_all(directory.join(format!("{name}_{ORTHODB}")).join(HMMS))?;
            }
        }
        Ok(Self {
            hmm,
            lineages,
            fingerprint: fold(&cutoffs),
        })
    }

    pub fn fingerprint(&self) -> u64 {
        self.fingerprint
    }

    fn bars(&self, lineage: &Lineage) -> Bars {
        lineage
            .models
            .iter()
            .map(|(model, cutoff)| {
                let bar = Bar {
                    sequence: *cutoff,
                    domain: *cutoff,
                };
                (model.clone(), bar)
            })
            .collect()
    }

    // Introns split one gene into several calls, and each call scores only its share of the cutoff.
    pub(super) fn search(
        &self,
        proteins: &Path,
        called: &[Orf],
        set: &MarkerSet,
        span: f64,
    ) -> Result<Vec<(usize, Hit)>> {
        let every = self
            .lineages
            .iter()
            .flat_map(|lineage| self.bars(lineage))
            .collect::<Bars>();
        let table = hmmer_engine::search(&self.hmm, proteins, fragments::floor(&every, span))?;
        let mut found = Vec::new();
        for lineage in &self.lineages {
            let bars = self.bars(lineage);
            let mut hits = fragments::complete(&table, &bars);
            let pieces =
                fragments::accepted(&table, &bars, span, |protein| !hits.contains_key(&protein));
            hits.extend(pieces);
            for (protein, best) in hits {
                let (Some(marker), Some(orf)) = (set.id(&best.model), called.get(protein)) else {
                    continue;
                };
                let mut hit = Hit::called(marker, orf, &best);
                hit.partial |= set.piece(marker, hit.place.score);
                found.push((orf.contig, hit));
            }
        }
        Ok(found)
    }
}

fn run(command: &mut Command) -> Result<()> {
    let status = command
        .status()
        .with_context(|| format!("could not run {command:?}"))?;
    if !status.success() {
        bail!("{command:?} failed with {status}");
    }
    Ok(())
}

fn fetch(directory: &Path, name: &str) -> Result<()> {
    let folder = format!("{name}_{ORTHODB}");
    let url = format!("{SOURCE}/{folder}.{RELEASE}.tar.gz");
    info!("Eukaryote panels: fetching {url}");
    let staging = tempfile::tempdir_in(directory)?;
    let archive = staging.path().join("lineage.tar.gz");
    run(Command::new("curl")
        .args(["-fsSL", "--retry", "3", "-o"])
        .arg(&archive)
        .arg(&url))?;
    run(Command::new("tar")
        .arg("-xzf")
        .arg(&archive)
        .arg("-C")
        .arg(staging.path())
        .arg(format!("{folder}/{HMMS}"))
        .arg(format!("{folder}/{CUTOFFS}"))
        .arg(format!("{folder}/dataset.cfg")))?;
    let target = directory.join(&folder);
    if target.exists() {
        std::fs::remove_dir_all(&target)?;
    }
    std::fs::rename(staging.path().join(&folder), &target)?;
    Ok(())
}

// BUSCO names a model by its group number. The same number recurs across lineages.
fn join(directory: &Path, lineages: &[Lineage], target: &Path) -> Result<()> {
    crate::report_sink::write_atomically(target, |file| {
        let mut sink = BufWriter::new(file);
        for lineage in lineages {
            let hmms = directory
                .join(format!("{}_{ORTHODB}", lineage.name))
                .join(HMMS);
            for (model, _) in &lineage.models {
                let text = std::fs::read_to_string(hmms.join(format!("{model}.hmm")))?;
                for line in text.lines() {
                    match line.starts_with("NAME ") {
                        true => writeln!(sink, "NAME  {model}")?,
                        false => writeln!(sink, "{line}")?,
                    }
                }
            }
        }
        sink.flush()?;
        Ok(())
    })
}
