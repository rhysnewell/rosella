use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use log::{debug, warn};

use crate::digest::fold;
use crate::markers::checkm::Copies;
use crate::markers::replicon::Shape;
use crate::markers::{Hit, MarkerSet};

const FORMAT: &str = "rosella-markers-7";

// Version 6 searched CheckM's models in the same pass, so its copies were settled on E-values of
// another search size. Its GTDB hits still hold and its CheckM copies read as never searched.
const SHARED_PASS_FORMAT: &str = "rosella-markers-6";

const UNSEARCHED: &str = "-";

// Bump when the annotation this file holds would come out different, whether that is what
// the search is handed or how a protein is settled between two models afterwards.
const FRAGMENT_PASS: u32 = 3;

const PATH_FIELD: usize = 3;

// An entry built at a floor holds every contig at or above it, so it serves any higher cutoff.
const FLOOR_FIELD: usize = 5;

// An entry from before the modified time joined the key still matches on the rest, since
// re-annotating every cached assembly costs more than the rare rewrite in place it would catch.
const STAMP_FIELD: usize = 9;

const ENTRY_PREFIX: &str = "markers.";
const ENTRY_SUFFIX: &str = ".tsv";

#[derive(Default)]
pub struct Rows {
    pub names: Vec<String>,
    pub lengths: Vec<usize>,
    pub hits: Vec<Vec<Hit>>,
    pub shapes: Vec<Shape>,
    pub checkm: Vec<Option<Vec<Copies>>>,
}

impl Rows {
    pub fn at_least(self, floor: usize) -> Self {
        let keep = self
            .lengths
            .iter()
            .map(|length| *length >= floor)
            .collect::<Vec<_>>();
        Self {
            names: kept(self.names, &keep),
            lengths: kept(self.lengths, &keep),
            hits: kept(self.hits, &keep),
            shapes: kept(self.shapes, &keep),
            checkm: kept(self.checkm, &keep),
        }
    }

    pub fn fill_from(mut self, held: Self, ceiling: usize) -> Result<Self> {
        let index = held
            .names
            .iter()
            .enumerate()
            .map(|(at, name)| (name.as_str(), at))
            .collect::<HashMap<_, _>>();
        let mut hits = held.hits;
        let mut checkm = held.checkm;
        for at in (0..self.names.len()).filter(|at| self.lengths[*at] >= ceiling) {
            let Some(&from) = index.get(self.names[at].as_str()) else {
                bail!("the cached annotation does not hold {}", self.names[at]);
            };
            self.hits[at] = std::mem::take(&mut hits[from]);
            self.shapes[at] = held.shapes[from];
            self.checkm[at] = std::mem::take(&mut checkm[from]);
        }
        Ok(self)
    }
}

fn kept<T>(values: Vec<T>, keep: &[bool]) -> Vec<T> {
    values
        .into_iter()
        .zip(keep)
        .filter_map(|(value, keep)| keep.then_some(value))
        .collect()
}

// The ingredients live in the file rather than only in its name, so changing how the key is
// spelled never discards an annotation that is still correct.
pub fn key(assembly: &str, floor: usize, fragment_span: f64) -> Result<String> {
    let source = fs::metadata(assembly)?;
    Ok([
        env!("ROSELLA_GENE_CALLER").to_string(),
        FRAGMENT_PASS.to_string(),
        format!(
            "{:016x}+{:016x}",
            fold(crate::markers::HMM_GZ),
            crate::markers::checkm::fingerprint()
        ),
        settled(assembly),
        source.len().to_string(),
        floor.to_string(),
        // The two gene rules were deleted at their defaults. Their zeros stay in the key so
        // annotations cached before that are still found.
        "0".to_string(),
        "0".to_string(),
        format!("{:016x}", fragment_span.to_bits()),
        stamp(&source),
    ]
    .join("\t"))
}

fn stamp(source: &fs::Metadata) -> String {
    source
        .modified()
        .ok()
        .and_then(|at| at.duration_since(std::time::UNIX_EPOCH).ok())
        .map_or(0, |since| since.as_nanos())
        .to_string()
}

fn settled(path: &str) -> String {
    fs::canonicalize(path)
        .map(|found| found.to_string_lossy().into_owned())
        .unwrap_or_else(|_| path.to_string())
}

// One assembly reaches rosella under many spellings, and re-annotating is most of a run, so
// the path is the one field compared through the filesystem rather than byte for byte.
fn held_floor(wanted: &str, held: &str) -> Option<usize> {
    let wanted = wanted.split('\t').collect::<Vec<_>>();
    let held = held.split('\t').collect::<Vec<_>>();
    if (held.len() != wanted.len() && held.len() != STAMP_FIELD)
        || wanted
            .iter()
            .zip(&held)
            .enumerate()
            .any(|(at, (ours, theirs))| at != PATH_FIELD && at != FLOOR_FIELD && ours != theirs)
    {
        return None;
    }
    (wanted[PATH_FIELD] == held[PATH_FIELD] || settled(held[PATH_FIELD]) == wanted[PATH_FIELD])
        .then(|| held[FLOOR_FIELD].parse().ok())
        .flatten()
}

fn named(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| name.starts_with(ENTRY_PREFIX) && name.ends_with(ENTRY_SUFFIX))
}

// The highest floor at or under the cutoff reads the fewest rows. Failing that, the lowest floor
// above it leaves the fewest contigs to call.
pub fn find(directory: &Path, key: &str) -> Option<(PathBuf, usize)> {
    let cutoff = key.split('\t').nth(FLOOR_FIELD)?.parse::<usize>().ok()?;
    let mut seen = 0usize;
    let mut ours = 0usize;
    let mut found: Option<(PathBuf, usize)> = None;
    for entry in fs::read_dir(directory).into_iter().flatten().flatten() {
        seen += 1;
        let path = entry.path();
        if !named(&path) {
            continue;
        }
        ours += 1;
        let Some(floor) = header(&path).and_then(|held| held_floor(key, &held)) else {
            continue;
        };
        let distance = |floor: usize| (floor > cutoff, floor.abs_diff(cutoff));
        if found
            .as_ref()
            .is_none_or(|(_, best)| distance(floor) < distance(*best))
        {
            found = Some((path, floor));
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
    let line = line.trim_end();
    line.strip_prefix(FORMAT)
        .or_else(|| line.strip_prefix(SHARED_PASS_FORMAT))?
        .strip_prefix('\t')
        .map(str::to_string)
}

pub fn record(
    directory: &Path,
    key: &str,
    set: &MarkerSet,
    found: &HashMap<&str, &[Copies]>,
) -> Result<()> {
    let Some((path, _)) = find(directory, key) else {
        return Ok(());
    };
    let Some(held) = header(&path) else {
        return Ok(());
    };
    let mut rows = read(&path, set)?;
    for (name, copies) in rows.names.iter().zip(&mut rows.checkm) {
        if let Some(searched) = found.get(name.as_str()) {
            *copies = Some(searched.to_vec());
        }
    }
    write(&path, &held, set, &rows)
}

pub fn read(path: &Path, set: &MarkerSet) -> Result<Rows> {
    let mut lines = BufReader::new(fs::File::open(path)?).lines();
    let Some(first) = lines.next().transpose()? else {
        bail!("{} is empty", path.display());
    };
    let shared_pass = first.starts_with(SHARED_PASS_FORMAT);
    if !first.starts_with(FORMAT) && !shared_pass {
        bail!("{} is not a marker cache", path.display());
    }
    let mut rows = Rows::default();
    for line in lines {
        let line = line?;
        let mut fields = line.splitn(6, '\t');
        let name = fields.next().unwrap_or_default();
        let hits = fields.next().unwrap_or_default();
        let coding_bases = fields.next().and_then(|field| field.parse().ok());
        let genes = fields.next().and_then(|field| field.parse().ok());
        let length = fields.next().and_then(|field| field.parse().ok());
        let checkm = fields.next().unwrap_or_default();
        let (Some(coding_bases), Some(genes), Some(length)) = (coding_bases, genes, length) else {
            bail!("{} has no gene shape or length for {name}", path.display());
        };
        rows.names.push(name.to_string());
        rows.lengths.push(length);
        rows.shapes.push(Shape {
            coding_bases,
            genes,
        });
        rows.hits.push(
            hits.split(',')
                .filter(|field| !field.is_empty())
                .filter_map(|field| Hit::decode(field, set))
                .collect(),
        );
        rows.checkm
            .push((!shared_pass && checkm != UNSEARCHED).then(|| {
                checkm
                    .split(',')
                    .filter_map(|field| decode_copies(field, set))
                    .collect()
            }));
    }
    Ok(rows)
}

fn decode_copies(field: &str, set: &MarkerSet) -> Option<Copies> {
    let mut parts = field.split(':');
    let lineage = set.checkm.lineage(parts.next()?)?;
    let model = set.checkm.id(parts.next()?)?;
    Some(Copies {
        set: u8::try_from(lineage).ok()?,
        model,
        copies: parts.next()?.parse().ok()?,
    })
}

pub fn write(path: &Path, key: &str, set: &MarkerSet, rows: &Rows) -> Result<()> {
    fs::create_dir_all(path.parent().unwrap_or(Path::new(".")))?;
    crate::report_sink::write_atomically(path, |file| entries(file, key, set, rows))
}

fn entries(file: &fs::File, key: &str, set: &MarkerSet, rows: &Rows) -> Result<()> {
    let mut sink = BufWriter::new(file);
    writeln!(sink, "{FORMAT}\t{key}")?;
    let each = rows
        .names
        .iter()
        .zip(&rows.hits)
        .zip(&rows.shapes)
        .zip(&rows.lengths)
        .zip(&rows.checkm);
    for ((((name, hits), shape), length), checkm) in each {
        write!(sink, "{name}\t")?;
        for (at, hit) in hits.iter().enumerate() {
            let separator = if at == 0 { "" } else { "," };
            write!(sink, "{separator}{}", hit.encode(set))?;
        }
        write!(
            sink,
            "\t{}\t{}\t{length}\t",
            shape.coding_bases, shape.genes
        )?;
        let Some(checkm) = checkm else {
            writeln!(sink, "{UNSEARCHED}")?;
            continue;
        };
        for (at, entry) in checkm.iter().enumerate() {
            let separator = if at == 0 { "" } else { "," };
            write!(
                sink,
                "{separator}{}:{}:{}",
                set.checkm.lineage_name(usize::from(entry.set)),
                set.checkm.name(entry.model),
                entry.copies
            )?;
        }
        writeln!(sink)?;
    }
    sink.flush()?;
    Ok(())
}
