use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::Result;

use super::Call;
use crate::digest::fold;
use crate::markers::Annotator;
use crate::markers::cache::{settled, stamp};

const FORMAT: &str = "rosella-elements-2";

// Calls join contigs by name, so a run over other bins searches only the contigs it adds.
pub fn key(annotator: &Annotator, models: u64) -> Result<String> {
    let assembly = match annotator.content {
        Some(content) => format!("{content:016x}"),
        None => {
            let source = std::fs::metadata(&annotator.assembly)?;
            format!(
                "{}\t{}\t{}",
                settled(&annotator.assembly),
                source.len(),
                stamp(&source)
            )
        }
    };
    Ok(format!(
        "{}\t{models:016x}\t{assembly}",
        env!("ROSELLA_GENE_CALLER")
    ))
}

fn path(directory: &Path, key: &str) -> PathBuf {
    directory.join(format!("elements.{:016x}.tsv", fold(key.as_bytes())))
}

pub fn read(directory: &Path, key: &str) -> HashMap<String, Call> {
    let Ok(file) = std::fs::File::open(path(directory, key)) else {
        return HashMap::new();
    };
    let mut lines = BufReader::new(file).lines().map_while(Result::ok);
    if lines.next() != Some(format!("{FORMAT}\t{key}")) {
        return HashMap::new();
    }
    lines
        .filter_map(|line| {
            let mut fields = line.split('\t');
            let name = fields.next()?.to_string();
            let call = Call {
                virus: fields.next()?.parse().ok()?,
                plasmid: fields.next()?.parse().ok()?,
            };
            Some((name, call))
        })
        .collect()
}

pub fn write(directory: &Path, key: &str, calls: &HashMap<String, Call>) -> Result<()> {
    std::fs::create_dir_all(directory)?;
    let mut names = calls.keys().collect::<Vec<_>>();
    names.sort_unstable();
    crate::report_sink::write_atomically(&path(directory, key), |file| {
        let mut sink = BufWriter::new(file);
        writeln!(sink, "{FORMAT}\t{key}")?;
        for name in names {
            let call = calls[name];
            writeln!(sink, "{name}\t{}\t{}", call.virus, call.plasmid)?;
        }
        sink.flush()?;
        Ok(())
    })
}
