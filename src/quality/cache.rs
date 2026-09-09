use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};

use super::tables::METADATA;
use crate::get_file_reader;

const FORMAT: &str = "rosella-gene-families-1";

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

/// The tables depend on the sequence, the contigs that survived the length filter and the
/// database, and on nothing else, so those three are the whole key.
pub fn path_for(directory: &Path, assembly: &str, names: &[String], database: &Path) -> PathBuf {
    let mut state = 0xcbf2_9ce4_8422_2325u64;
    mix(&mut state, assembly.as_bytes());
    mix(&mut state, database.to_string_lossy().as_bytes());
    if let Ok(stats) = std::fs::metadata(assembly) {
        mix(&mut state, &stats.len().to_le_bytes());
    }
    mix(&mut state, &names.len().to_le_bytes());
    for name in names {
        mix(&mut state, name.as_bytes());
        mix(&mut state, b"\0");
    }
    directory.join(format!("{state:016x}.families.gz"))
}

pub fn read(path: &Path, contigs: usize) -> Result<Annotation> {
    let mut reader = flate2::read::GzDecoder::new(get_file_reader(path)?);
    let mut text = String::new();
    std::io::Read::read_to_string(&mut reader, &mut text)?;

    let mut lines = text.lines();
    match lines.next() {
        Some(FORMAT) => {}
        _ => bail!("{} was not written by this version", path.display()),
    }

    let mut metadata = Vec::with_capacity(contigs);
    let mut hits = Vec::with_capacity(contigs);
    for line in lines {
        let (counts, families) = line.split_once('\t').unwrap_or((line, ""));
        let mut row = [0u32; METADATA];
        let mut columns = counts.split(' ');
        for cell in row.iter_mut() {
            *cell = columns
                .next()
                .and_then(|value| value.parse().ok())
                .ok_or_else(|| anyhow!("{} has a short metadata row", path.display()))?;
        }
        let mut held = Vec::new();
        for pair in families.split(' ').filter(|pair| !pair.is_empty()) {
            let (family, count) = pair
                .split_once(':')
                .ok_or_else(|| anyhow!("{} has an unreadable gene family", path.display()))?;
            held.push((family.parse()?, count.parse()?));
        }
        metadata.push(row);
        hits.push(held);
    }

    if metadata.len() != contigs {
        bail!(
            "{} holds {} contigs against the {} handed in",
            path.display(),
            metadata.len(),
            contigs
        );
    }
    Ok(Annotation { metadata, hits })
}

pub fn write(path: &Path, annotation: &Annotation) -> Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut sink = BufWriter::new(flate2::write::GzEncoder::new(
        std::fs::File::create(path)?,
        flate2::Compression::default(),
    ));
    writeln!(sink, "{FORMAT}")?;
    for (row, held) in annotation.metadata.iter().zip(&annotation.hits) {
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
