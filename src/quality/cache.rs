use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};

use super::tables::METADATA;
use crate::get_file_reader;

const FORMAT: &str = "rosella-gene-families-2";

pub struct Annotation {
    pub metadata: Vec<[u32; METADATA]>,
    pub hits: Vec<Vec<(u32, u32)>>,
}

fn mix(state: &mut u64, bytes: &[u8]) {
    for byte in bytes {
        *state ^= *byte as u64;
        *state = state.wrapping_mul(0x0000_0100_0000_01b3);
    }
}

/// The tables depend on the sequence and the database and on nothing else. Rows carry their
/// contig name, so one file serves any contig length floor over the same assembly.
pub fn path_for(directory: &Path, assembly: &str, database: &Path) -> PathBuf {
    let mut state = 0xcbf2_9ce4_8422_2325u64;
    mix(&mut state, assembly.as_bytes());
    mix(&mut state, database.to_string_lossy().as_bytes());
    if let Ok(stats) = std::fs::metadata(assembly) {
        mix(&mut state, &stats.len().to_le_bytes());
    }
    directory.join(format!("{state:016x}.families.gz"))
}

type Row = (String, [u32; METADATA], Vec<(u32, u32)>);

fn row(line: &str, path: &Path) -> Result<Row> {
    let mut fields = line.split('\t');
    let Some(name) = fields.next() else {
        bail!("{} has a row with no contig", path.display());
    };
    let mut counts = fields.next().unwrap_or_default().split(' ');
    let mut metadata = [0u32; METADATA];
    for cell in metadata.iter_mut() {
        *cell = counts
            .next()
            .and_then(|value| value.parse().ok())
            .ok_or_else(|| anyhow!("{} has a short metadata row", path.display()))?;
    }
    let mut held = Vec::new();
    for pair in fields
        .next()
        .unwrap_or_default()
        .split(' ')
        .filter(|pair| !pair.is_empty())
    {
        let (family, count) = pair
            .split_once(':')
            .ok_or_else(|| anyhow!("{} has an unreadable gene family", path.display()))?;
        held.push((family.parse()?, count.parse()?));
    }
    Ok((name.to_string(), metadata, held))
}

pub fn read(path: &Path) -> Result<(Vec<String>, Annotation)> {
    let mut reader = flate2::read::GzDecoder::new(get_file_reader(path)?);
    let mut text = String::new();
    std::io::Read::read_to_string(&mut reader, &mut text)?;

    let mut lines = text.lines();
    match lines.next() {
        Some(FORMAT) => {}
        _ => bail!("{} was not written by this version", path.display()),
    }

    let mut names = Vec::new();
    let mut metadata = Vec::new();
    let mut hits = Vec::new();
    for line in lines {
        let (name, row, found) = row(line, path)?;
        names.push(name);
        metadata.push(row);
        hits.push(found);
    }
    Ok((names, Annotation { metadata, hits }))
}

pub fn select(names: &[String], held: &[String], annotation: Annotation) -> Result<Annotation> {
    let index = held
        .iter()
        .enumerate()
        .map(|(position, name)| (name.as_str(), position))
        .collect::<HashMap<_, _>>();
    let mut metadata = Vec::with_capacity(names.len());
    let mut hits = Vec::with_capacity(names.len());
    let Annotation {
        metadata: rows,
        hits: found,
    } = annotation;
    for name in names {
        let Some(position) = index.get(name.as_str()) else {
            bail!("the gene family tables do not hold {name}");
        };
        metadata.push(rows[*position]);
        hits.push(found[*position].clone());
    }
    Ok(Annotation { metadata, hits })
}

pub fn write(path: &Path, names: &[String], annotation: &Annotation) -> Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut sink = BufWriter::new(flate2::write::GzEncoder::new(
        std::fs::File::create(path)?,
        flate2::Compression::default(),
    ));
    writeln!(sink, "{FORMAT}")?;
    for ((name, row), held) in names.iter().zip(&annotation.metadata).zip(&annotation.hits) {
        write!(sink, "{name}\t")?;
        for (column, count) in row.iter().enumerate() {
            match column {
                0 => write!(sink, "{count}")?,
                _ => write!(sink, " {count}")?,
            }
        }
        sink.write_all(b"\t")?;
        for (position, (family, count)) in held.iter().enumerate() {
            match position {
                0 => write!(sink, "{family}:{count}")?,
                _ => write!(sink, " {family}:{count}")?,
            }
        }
        sink.write_all(b"\n")?;
    }
    sink.flush()?;
    Ok(())
}
