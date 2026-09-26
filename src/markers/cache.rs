use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use log::{debug, warn};

use crate::digest::fold;
use crate::markers::replicon::Shape;
use crate::markers::{Hit, MarkerSet};

const FORMAT: &str = "rosella-markers-4";

/// Bump when the annotation this file holds would come out different, whether that is what
/// the search is handed or how a protein is settled between two models afterwards.
const FRAGMENT_PASS: u32 = 3;

const PATH_FIELD: usize = 3;

// An entry built at a floor holds every contig at or above it, so it serves any higher cutoff.
const FLOOR_FIELD: usize = 5;

const ENTRY_PREFIX: &str = "markers.";
const ENTRY_SUFFIX: &str = ".tsv";

#[derive(Default)]
pub struct Rows {
    pub names: Vec<String>,
    pub lengths: Vec<usize>,
    pub hits: Vec<Vec<Hit>>,
    pub shapes: Vec<Shape>,
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
        for at in (0..self.names.len()).filter(|at| self.lengths[*at] >= ceiling) {
            let Some(&from) = index.get(self.names[at].as_str()) else {
                bail!("the cached annotation does not hold {}", self.names[at]);
            };
            self.hits[at] = std::mem::take(&mut hits[from]);
            self.shapes[at] = held.shapes[from];
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

/// The ingredients live in the file rather than only in its name, so changing how the key is
/// spelled never discards an annotation that is still correct.
pub fn key(assembly: &str, floor: usize, fragment_span: f64) -> Result<String> {
    let source = fs::metadata(assembly)?;
    Ok([
        env!("ROSELLA_GENE_CALLER").to_string(),
        FRAGMENT_PASS.to_string(),
        format!("{:016x}", fold(crate::markers::HMM_GZ)),
        settled(assembly),
        source.len().to_string(),
        floor.to_string(),
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
fn held_floor(wanted: &str, held: &str) -> Option<usize> {
    let wanted = wanted.split('\t').collect::<Vec<_>>();
    let held = held.split('\t').collect::<Vec<_>>();
    if wanted.len() != held.len()
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
    line.trim_end()
        .strip_prefix(FORMAT)?
        .strip_prefix('\t')
        .map(str::to_string)
}

pub fn read(path: &Path, set: &MarkerSet) -> Result<Rows> {
    let mut lines = BufReader::new(fs::File::open(path)?).lines();
    let Some(first) = lines.next().transpose()? else {
        bail!("{} is empty", path.display());
    };
    if !first.starts_with(FORMAT) {
        bail!("{} is not a marker cache", path.display());
    }
    let mut rows = Rows::default();
    for line in lines {
        let line = line?;
        let mut fields = line.splitn(5, '\t');
        let name = fields.next().unwrap_or_default();
        let hits = fields.next().unwrap_or_default();
        let coding_bases = fields.next().and_then(|field| field.parse().ok());
        let genes = fields.next().and_then(|field| field.parse().ok());
        let length = fields.next().and_then(|field| field.parse().ok());
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
    Ok(rows)
}

pub fn write(path: &Path, key: &str, set: &MarkerSet, rows: &Rows) -> Result<()> {
    let parent = path.parent().unwrap_or(Path::new("."));
    fs::create_dir_all(parent)?;
    // A uniquely named file beside the entry, renamed into place, so two processes annotating
    // one assembly never interleave writes and the last complete file wins.
    let pending = tempfile::NamedTempFile::new_in(parent)?;
    let mut sink = BufWriter::new(pending.as_file());
    writeln!(sink, "{FORMAT}\t{key}")?;
    let each = rows
        .names
        .iter()
        .zip(&rows.hits)
        .zip(&rows.shapes)
        .zip(&rows.lengths);
    for (((name, hits), shape), length) in each {
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
        writeln!(sink, "\t{}\t{}\t{length}", shape.coding_bases, shape.genes)?;
    }
    sink.flush()?;
    drop(sink);
    pending.persist(path)?;
    Ok(())
}
