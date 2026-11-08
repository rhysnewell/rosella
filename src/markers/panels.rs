use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

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
const CUTOFFS: &str = "scores_cutoff";
const HMMS: &str = "hmms";

pub struct Lineage {
    pub name: &'static str,
    pub models: Vec<(String, f64)>,
    hmm: PathBuf,
}

impl Lineage {
    fn bars(&self) -> Bars {
        self.models
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
}

pub struct Panel {
    pub lineages: Vec<Lineage>,
    fingerprint: u64,
}

impl Panel {
    pub fn prepare(directory: &Path, names: &[&'static str], fetch: bool) -> Result<Self> {
        std::fs::create_dir_all(directory)?;
        let mut lineages = Vec::new();
        let mut cutoffs = Vec::new();
        for &name in names {
            let folder = directory.join(format!("{name}_{ORTHODB}"));
            let hmm = directory.join(format!("{name}.hmm"));
            if !folder.join(CUTOFFS).exists() || !(hmm.exists() || folder.join(HMMS).exists()) {
                if !fetch {
                    bail!(
                        "the BUSCO {name} set is not downloaded, so run once with --download-data <DIR>"
                    );
                }
                download(directory, name)?;
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
            let lineage = Lineage { name, models, hmm };
            if !lineage.hmm.exists() {
                join(&folder.join(HMMS), &lineage)?;
                std::fs::remove_dir_all(folder.join(HMMS))?;
            }
            lineages.push(lineage);
        }
        Ok(Self {
            lineages,
            fingerprint: fold(&cutoffs),
        })
    }

    pub fn fingerprint(&self) -> u64 {
        self.fingerprint
    }

    // Introns split one gene into several calls, and each call scores only its share of the cutoff.
    pub(super) fn search(
        &self,
        proteins: &Path,
        called: &[Orf],
        set: &MarkerSet,
        span: f64,
    ) -> Result<Vec<(usize, Hit)>> {
        let mut found = Vec::new();
        for lineage in &self.lineages {
            let bars = lineage.bars();
            let table =
                hmmer_engine::search(&lineage.hmm, proteins, fragments::floor(&bars, span))?;
            let (hits, _) = fragments::called(&table, &bars, span, |_| true);
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

fn download(directory: &Path, name: &str) -> Result<()> {
    let folder = format!("{name}_{ORTHODB}");
    let url = format!("{SOURCE}/{folder}.{RELEASE}.tar.gz");
    info!("Eukaryote panels: fetching {url}");
    let staging = tempfile::tempdir_in(directory)?;
    let body = ureq::get(&url)
        .call()
        .with_context(|| format!("could not fetch {url}"))?
        .into_body();
    let hmms = Path::new(&folder).join(HMMS);
    let kept = [CUTOFFS, "dataset.cfg"].map(|file| Path::new(&folder).join(file));
    let mut archive = tar::Archive::new(flate2::read::GzDecoder::new(body.into_reader()));
    for entry in archive.entries()? {
        let mut entry = entry?;
        let path = entry.path()?.into_owned();
        if path.starts_with(&hmms) || kept.contains(&path) {
            entry.unpack_in(staging.path())?;
        }
    }
    let target = directory.join(&folder);
    if target.exists() {
        std::fs::remove_dir_all(&target)?;
    }
    std::fs::rename(staging.path().join(&folder), &target)?;
    Ok(())
}

// A BUSCO model file names the model by its bare group number, the cutoff table by its full id.
fn join(hmms: &Path, lineage: &Lineage) -> Result<()> {
    crate::report_sink::write_atomically(&lineage.hmm, |file| {
        let mut sink = BufWriter::new(file);
        for (model, _) in &lineage.models {
            let text = std::fs::read_to_string(hmms.join(format!("{model}.hmm")))?;
            for line in text.lines() {
                match line.starts_with("NAME ") {
                    true => writeln!(sink, "NAME  {model}")?,
                    false => writeln!(sink, "{line}")?,
                }
            }
        }
        sink.flush()?;
        Ok(())
    })
}
