use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use log::{debug, warn};

use crate::markers::{Hit, MarkerSet};

const FORMAT: &str = "rosella-markers-2";

/// Bump when the fragment rescue or the protein filter changes what the search is handed.
const FRAGMENT_PASS: u32 = 3;

/// Stable across compiler versions, unlike the hasher in std, so a rebuild does not rename
/// every entry.
fn fold(bytes: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for byte in bytes {
        hash ^= *byte as u64;
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

const PATH_FIELD: usize = 3;

const ENTRY_PREFIX: &str = "markers.";
const ENTRY_SUFFIX: &str = ".tsv";

/// The ingredients live in the file rather than only in its name, so changing how the key is
/// spelled never discards an annotation that is still correct.
pub fn key(assembly: &str, min_contig_size: usize, fragment_span: f64) -> Result<String> {
    let source = fs::metadata(assembly)?;
    Ok([
        env!("ROSELLA_GENE_CALLER").to_string(),
        FRAGMENT_PASS.to_string(),
        format!("{:016x}", fold(crate::markers::HMM_GZ)),
        settled(assembly),
        source.len().to_string(),
        min_contig_size.to_string(),
        // The two gene rules were deleted at their defaults. Their zeros stay in the key so
        // annotations cached before that are still found.
        "0".to_string(),
        "0".to_string(),
        format!("{:016x}", fragment_span.to_bits()),
    ]
    .join("\t"))
}

fn settled(path: &str) -> String {
    fs::canonicalize(path)
        .map(|found| found.to_string_lossy().into_owned())
        .unwrap_or_else(|_| path.to_string())
}

/// One assembly reaches rosella under many spellings, and re-annotating is most of a run, so
/// the path is the one field compared through the filesystem rather than byte for byte.
fn same(wanted: &str, held: &str) -> bool {
    let wanted = wanted.split('\t').collect::<Vec<_>>();
    let held = held.split('\t').collect::<Vec<_>>();
    if wanted.len() != held.len() {
        return false;
    }
    if wanted
        .iter()
        .zip(&held)
        .enumerate()
        .any(|(at, (ours, theirs))| at != PATH_FIELD && ours != theirs)
    {
        return false;
    }
    wanted[PATH_FIELD] == held[PATH_FIELD] || settled(held[PATH_FIELD]) == wanted[PATH_FIELD]
}

fn named(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| name.starts_with(ENTRY_PREFIX) && name.ends_with(ENTRY_SUFFIX))
}

pub fn find(directory: &Path, key: &str) -> Option<PathBuf> {
    let mut seen = 0usize;
    let mut ours = 0usize;
    let mut found = None;
    for entry in fs::read_dir(directory).into_iter().flatten().flatten() {
        seen += 1;
        let path = entry.path();
        if !named(&path) {
            continue;
        }
        ours += 1;
        if found.is_none() && header(&path).is_some_and(|held| same(key, &held)) {
            found = Some(path);
        }
    }
    if found.is_none() {
        if seen > 0 && ours == 0 {
            warn!(
                "{} holds no marker cache entry. Is it the marker database rather than the cache?",
                directory.display()
            );
        }
        debug!(
            "No marker annotation cached in {} for {key}",
            directory.display()
        );
    }
    found
}

pub fn write_path(directory: &Path, key: &str) -> PathBuf {
    directory.join(format!(
        "{ENTRY_PREFIX}{:016x}{ENTRY_SUFFIX}",
        fold(key.as_bytes())
    ))
}

fn header(path: &Path) -> Option<String> {
    let mut line = String::new();
    BufReader::new(fs::File::open(path).ok()?)
        .read_line(&mut line)
        .ok()?;
    line.trim_end()
        .strip_prefix(FORMAT)?
        .strip_prefix('\t')
        .map(str::to_string)
}

pub fn read(path: &Path, set: &MarkerSet) -> Result<(Vec<String>, Vec<Vec<Hit>>)> {
    let mut lines = BufReader::new(fs::File::open(path)?).lines();
    let Some(first) = lines.next().transpose()? else {
        bail!("{} is empty", path.display());
    };
    if !first.starts_with(FORMAT) {
        bail!("{} is not a marker cache", path.display());
    }
    let mut names = Vec::new();
    let mut per_contig = Vec::new();
    for line in lines {
        let line = line?;
        let (name, hits) = line.split_once('\t').unwrap_or((line.as_str(), ""));
        names.push(name.to_string());
        per_contig.push(
            hits.split(',')
                .filter(|field| !field.is_empty())
                .filter_map(|field| {
                    let (model, partial) = field.split_once(':')?;
                    Some(Hit {
                        marker: set.id(model)?,
                        partial: partial == "1",
                    })
                })
                .collect(),
        );
    }
    Ok((names, per_contig))
}

pub fn write(
    path: &Path,
    key: &str,
    set: &MarkerSet,
    names: &[String],
    per_contig: &[Vec<Hit>],
) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    // Written beside the final name and renamed, so two arms racing one cache cannot leave a
    // half file that the next run would read as the whole annotation.
    let pending = path.with_extension("pending");
    let mut sink = BufWriter::new(fs::File::create(&pending)?);
    writeln!(sink, "{FORMAT}\t{key}")?;
    for (name, hits) in names.iter().zip(per_contig) {
        write!(sink, "{name}\t")?;
        for (at, hit) in hits.iter().enumerate() {
            let separator = if at == 0 { "" } else { "," };
            write!(
                sink,
                "{separator}{}:{}",
                set.name(hit.marker),
                u8::from(hit.partial)
            )?;
        }
        writeln!(sink)?;
    }
    sink.flush()?;
    drop(sink);
    fs::rename(&pending, path)?;
    Ok(())
}
