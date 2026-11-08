use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use log::debug;

use super::Families;
use crate::external::hmmer_engine;
use crate::markers::fragments::cutoffs;

pub const RELEASE: &str = "Pfam39.0";
const SOURCE: &str = "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases";
pub const MODELS: &str = "elements.hmm";
const BARS: &str = "gathering.tsv";

pub fn complete(folder: &Path) -> bool {
    [MODELS, BARS].iter().all(|file| folder.join(file).exists())
}

pub fn fetch(directory: &Path, families: &Families) -> Result<()> {
    std::fs::create_dir_all(directory)?;
    let url = format!("{SOURCE}/{RELEASE}/Pfam-A.hmm.gz");
    debug!("Element families: fetching {url}");
    let staging = tempfile::tempdir_in(directory)?;
    let folder = staging.path().join(RELEASE);
    std::fs::create_dir(&folder)?;
    let body = ureq::get(&url)
        .call()
        .with_context(|| format!("could not fetch {url}"))?
        .into_body();
    let reader = BufReader::new(flate2::read::MultiGzDecoder::new(body.into_reader()));
    let mut models = BufWriter::new(File::create(folder.join(MODELS))?);
    let mut bars = String::new();
    hmmer_engine::each_model(reader, |name, block| {
        let gathering = block
            .lines()
            .find_map(|line| line.strip_prefix("GA "))
            .and_then(cutoffs);
        let (true, Some((sequence, domain))) = (families.holds(name), gathering) else {
            return Ok(());
        };
        models.write_all(block.as_bytes())?;
        bars.push_str(&format!("{name}\t{sequence}\t{domain}\n"));
        Ok(())
    })?;
    models.flush()?;
    if bars.is_empty() {
        bail!("{url} held none of the element families");
    }
    debug!("Element families: kept {}", bars.lines().count());
    std::fs::write(folder.join(BARS), bars)?;
    let target = directory.join(RELEASE);
    match std::fs::rename(&folder, &target) {
        Ok(()) => Ok(()),
        Err(_) if complete(&target) => Ok(()),
        Err(error) => Err(error).with_context(|| format!("could not fill {}", target.display())),
    }
}

pub fn read_bars(folder: &Path) -> Result<HashMap<String, (f64, f64)>> {
    let path = folder.join(BARS);
    let text = std::fs::read_to_string(&path)?;
    text.lines()
        .map(|line| {
            let mut fields = line.split('\t');
            let (Some(name), Some(sequence), Some(domain)) =
                (fields.next(), fields.next(), fields.next())
            else {
                bail!("{} holds a bad row {line}", path.display());
            };
            Ok((name.to_string(), (sequence.parse()?, domain.parse()?)))
        })
        .collect()
}
