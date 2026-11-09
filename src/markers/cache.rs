use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use log::{debug, warn};

use crate::digest::fold;
use crate::markers::checkm::Copies;
use crate::markers::replicon::Shape;
use crate::markers::{Hit, MarkerSet};

// Version 10 marks each CheckM copy whose gene a contig end cuts, so no older entry holds what a
// run reads. Version 9 was a rolled back layout whose entries still sit in some caches.
const FORMAT: &str = "rosella-markers-10";

// Bump when the annotation this file holds would come out different, whether that is what
// the search is handed or how a protein is settled between two models afterwards.
const FRAGMENT_PASS: u32 = 3;

const PATH_FIELD: usize = 3;
const SIZE_FIELD: usize = 4;

// An entry built at a floor holds every contig at or above it, so it serves any higher cutoff.
const FLOOR_FIELD: usize = 5;

// An entry from before the modified time joined the key still matches on the rest, since
// re-annotating every cached assembly costs more than the rare rewrite in place it would catch.
const STAMP_FIELD: usize = 9;

// Rows are joined to contigs by name, so a copied, renamed or reordered assembly holding the
// same contigs reuses its entry. An entry without it falls back to the file's path, size and stamp.
const CONTENT_FIELD: usize = 10;

const UNKNOWN: &str = "-";

const FILE_FIELDS: [usize; 4] = [PATH_FIELD, SIZE_FIELD, STAMP_FIELD, CONTENT_FIELD];

const ENTRY_PREFIX: &str = "markers.";
const ENTRY_SUFFIX: &str = ".tsv";

#[derive(Default)]
pub struct Rows {
    pub names: Vec<String>,
    pub lengths: Vec<usize>,
    pub hits: Vec<Vec<Hit>>,
    pub shapes: Vec<Shape>,
    pub checkm: Vec<Vec<Copies>>,
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
        let index = crate::rows::positions(held.names.iter().map(String::as_str));
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

// The ingredients live in the file, so respelling the key never drops a good annotation. The
// panels join the model fingerprint, so an entry searched without them never serves one with them.
pub fn key(
    assembly: &str,
    content: Option<u64>,
    floor: usize,
    fragment_span: f64,
    panel: Option<u64>,
) -> Result<String> {
    let source = fs::metadata(assembly)?;
    let panel = panel.map_or_else(String::new, |panel| format!("+{panel:016x}"));
    let models = format!(
        "{:016x}+{:016x}{panel}",
        fold(crate::markers::HMM_GZ),
        crate::markers::checkm::fingerprint()
    );
    Ok([
        env!("ROSELLA_GENE_CALLER").to_string(),
        FRAGMENT_PASS.to_string(),
        models,
        settled(assembly),
        source.len().to_string(),
        floor.to_string(),
        // The two gene rules were deleted at their defaults. Their zeros stay in the key so
        // annotations cached before that are still found.
        "0".to_string(),
        "0".to_string(),
        format!("{:016x}", fragment_span.to_bits()),
        stamp(&source),
        content.map_or_else(|| UNKNOWN.to_string(), |content| format!("{content:016x}")),
    ]
    .join("\t"))
}

pub(crate) fn stamp(source: &fs::Metadata) -> String {
    source
        .modified()
        .ok()
        .and_then(|at| at.duration_since(std::time::UNIX_EPOCH).ok())
        .map_or(0, |since| since.as_nanos())
        .to_string()
}

pub(crate) fn settled(path: &str) -> String {
    fs::canonicalize(path)
        .map(|found| found.to_string_lossy().into_owned())
        .unwrap_or_else(|_| path.to_string())
}

// One assembly reaches rosella under many spellings, and re-annotating is most of a run, so
// the path is the one field compared through the filesystem rather than byte for byte.
fn held_floor(wanted: &str, held: &str) -> Option<(usize, bool)> {
    let wanted = wanted.split('\t').collect::<Vec<_>>();
    let held = held.split('\t').collect::<Vec<_>>();
    if ![STAMP_FIELD, CONTENT_FIELD, wanted.len()].contains(&held.len())
        || wanted
            .iter()
            .zip(&held)
            .enumerate()
            .any(|(at, (ours, theirs))| {
                at != FLOOR_FIELD && !FILE_FIELDS.contains(&at) && ours != theirs
            })
    {
        return None;
    }
    let same = match (content(&wanted), content(&held)) {
        (Some(ours), Some(theirs)) => ours == theirs,
        _ => same_file(&wanted, &held),
    };
    let floor = held[FLOOR_FIELD].parse().ok().filter(|_| same)?;
    Some((floor, content(&held).is_some()))
}

fn keyed(key: &str) -> bool {
    content(&key.split('\t').collect::<Vec<_>>()).is_some()
}

fn content<'a>(fields: &[&'a str]) -> Option<&'a str> {
    fields
        .get(CONTENT_FIELD)
        .copied()
        .filter(|content| *content != UNKNOWN)
}

fn same_file(wanted: &[&str], held: &[&str]) -> bool {
    (wanted[PATH_FIELD] == held[PATH_FIELD] || settled(held[PATH_FIELD]) == wanted[PATH_FIELD])
        && wanted[SIZE_FIELD] == held[SIZE_FIELD]
        && held
            .get(STAMP_FIELD)
            .is_none_or(|stamp| wanted.get(STAMP_FIELD) == Some(stamp))
}

// The entry is rewritten under the content key, so a copy of the assembly on another path or
// another machine finds it too.
pub fn upgrade(
    directory: &Path,
    key: &str,
    held: &Path,
    floor: usize,
    set: &MarkerSet,
    rows: &Rows,
) {
    if !keyed(key) || header(held).is_none_or(|header| keyed(&header)) {
        return;
    }
    let mut fields = key.split('\t').collect::<Vec<_>>();
    let floor = floor.to_string();
    fields[FLOOR_FIELD] = &floor;
    let key = fields.join("\t");
    let path = write_path(directory, &key);
    match write(&path, &key, set, rows) {
        Ok(()) => {
            debug!("Marker cache: rekeyed {}", path.display());
            if path != held {
                let _ = fs::remove_file(held);
            }
        }
        Err(error) => warn!("Marker cache: {error} ({})", path.display()),
    }
}

fn named(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| name.starts_with(ENTRY_PREFIX) && name.ends_with(ENTRY_SUFFIX))
}

// The highest floor at or under the cutoff reads the fewest rows, else the lowest above it calls
// the fewest contigs. At one floor a content key wins, so a path-keyed twin is never rekeyed again.
pub fn find(directory: &Path, key: &str) -> Option<(PathBuf, usize)> {
    let cutoff = key.split('\t').nth(FLOOR_FIELD)?.parse::<usize>().ok()?;
    let mut seen = 0usize;
    let mut ours = 0usize;
    let mut found: Option<(PathBuf, (usize, bool))> = None;
    for entry in fs::read_dir(directory).into_iter().flatten().flatten() {
        seen += 1;
        let path = entry.path();
        if !named(&path) {
            continue;
        }
        ours += 1;
        let Some(held) = header(&path).and_then(|held| held_floor(key, &held)) else {
            continue;
        };
        let distance =
            |(floor, keyed): (usize, bool)| (floor > cutoff, floor.abs_diff(cutoff), !keyed);
        if found
            .as_ref()
            .is_none_or(|(_, best)| distance(held) < distance(*best))
        {
            found = Some((path, held));
        }
    }
    if found.is_none() {
        static SAID: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
        if seen > 0 && ours == 0 && crate::progress::first(&SAID) {
            warn!(
                "Marker cache: no entry in {} (is it the marker database?)",
                directory.display()
            );
        }
        debug!(
            "Marker cache: no entry in {} for {key}",
            directory.display()
        );
    }
    found.map(|(path, (floor, _))| (path, floor))
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
    line.strip_prefix(FORMAT)?
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
        rows.checkm.push(
            checkm
                .split(',')
                .filter_map(|field| decode_copies(field, set))
                .collect(),
        );
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
        score: parts.next()?.parse().ok()?,
        partial: parts.next()? == "1",
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
        for (at, entry) in checkm.iter().enumerate() {
            let separator = if at == 0 { "" } else { "," };
            write!(
                sink,
                "{separator}{}:{}:{}:{}",
                set.checkm.lineage_name(usize::from(entry.set)),
                set.checkm.name(entry.model),
                entry.score,
                u8::from(entry.partial)
            )?;
        }
        writeln!(sink)?;
    }
    sink.flush()?;
    Ok(())
}
