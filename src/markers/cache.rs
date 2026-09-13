use std::fs;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};

use crate::markers::{Hit, MarkerSet};
use crate::quality::orfs::GeneRules;

const FORMAT: &str = "rosella-markers-1";

/// The build commit is part of the key, so a rebuilt binary never reads an annotation its own
/// gene caller or fragment pass would not produce.
pub fn path(
    directory: &Path,
    assembly: &str,
    min_contig_size: usize,
    genes: GeneRules,
    fragment_span: f64,
) -> Result<PathBuf> {
    let source = fs::metadata(assembly)?;
    let mut hasher = DefaultHasher::new();
    FORMAT.hash(&mut hasher);
    env!("ROSELLA_BUILD_COMMIT").hash(&mut hasher);
    assembly.hash(&mut hasher);
    source.len().hash(&mut hasher);
    min_contig_size.hash(&mut hasher);
    genes.min_length.hash(&mut hasher);
    genes.model_depth.hash(&mut hasher);
    fragment_span.to_bits().hash(&mut hasher);
    Ok(directory.join(format!("markers.{:016x}.tsv", hasher.finish())))
}

pub fn read(path: &Path, set: &MarkerSet) -> Result<(Vec<String>, Vec<Vec<Hit>>)> {
    let mut lines = BufReader::new(fs::File::open(path)?).lines();
    match lines.next().transpose()?.as_deref() {
        Some(FORMAT) => {}
        other => bail!("{} is not a marker cache: {other:?}", path.display()),
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
    writeln!(sink, "{FORMAT}")?;
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
