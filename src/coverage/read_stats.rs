use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use anyhow::Result;
use log::info;
use noodles_bam as bam;
use noodles_sam::alignment::record::{
    Cigar as _,
    data::field::{Tag, Value},
};

use super::depth_file::{Contigs, packed};

pub(super) struct ReadStats<'a> {
    names: &'a HashMap<&'a str, u32>,
    pub(super) reads: Vec<u32>,
    seen: u64,
    mapped: u64,
    tied: u64,
    read_bases: u64,
    mapped_bases: u64,
    links: HashMap<(u32, u32), u32>,
    ties: HashMap<(u32, u32), u32>,
}

impl<'a> ReadStats<'a> {
    pub(super) fn new(contigs: &'a Contigs, names: &'a HashMap<&'a str, u32>) -> Self {
        Self {
            names,
            reads: vec![0; contigs.names.len()],
            seen: 0,
            mapped: 0,
            tied: 0,
            read_bases: 0,
            mapped_bases: 0,
            links: HashMap::new(),
            ties: HashMap::new(),
        }
    }

    pub(super) fn note(&mut self, record: &bam::Record) -> Result<()> {
        let flags = record.flags();
        if flags.is_secondary() || flags.is_supplementary() {
            return Ok(());
        }
        self.seen += 1;
        let Some(contig) = record
            .reference_sequence_id()
            .transpose()?
            .filter(|_| !flags.is_unmapped())
        else {
            self.read_bases += record.sequence().len() as u64;
            return Ok(());
        };
        let cigar = record.cigar();
        self.mapped += 1;
        self.reads[contig] += 1;
        self.read_bases += cigar.read_length()? as u64;
        self.mapped_bases += cigar.alignment_span()? as u64;
        if let Some(Value::String(others)) =
            record.data().get(&Tag::OTHER_ALIGNMENTS).transpose()?
        {
            for other in others.split(|byte| *byte == b';') {
                let name = other.split(|byte| *byte == b',').next().unwrap_or_default();
                let Some(&linked) = std::str::from_utf8(name)
                    .ok()
                    .and_then(|name| self.names.get(name))
                else {
                    continue;
                };
                if linked as usize != contig {
                    *self
                        .links
                        .entry(ordered(contig as u32, linked))
                        .or_default() += 1;
                }
            }
        }
        Ok(())
    }

    pub(super) fn tie(&mut self, contigs: impl Iterator<Item = u32> + Clone) {
        self.tied += 1;
        for (at, first) in contigs.clone().enumerate() {
            for second in contigs.clone().skip(at + 1) {
                *self.ties.entry(ordered(first, second)).or_default() += 1;
            }
        }
    }

    pub(super) fn write(&self, directory: &Path, name: &str) -> Result<()> {
        let share = |part: u64, whole: u64| part as f64 / whole.max(1) as f64;
        let length = share(self.read_bases, self.seen);
        info!(
            "{name}: {} reads, {:.1}% mapped, {:.1}% tied, mean length {length:.0} bp",
            self.seen,
            100.0 * share(self.mapped, self.seen),
            100.0 * share(self.tied, self.mapped),
        );
        let summary = format!(
            "sample\treads\tmapped\ttied\tmean_read_length\tmapped_bases\n{name}\t{}\t{}\t{}\t{length:.1}\t{}\n",
            self.seen, self.mapped, self.tied, self.mapped_bases
        );
        crate::report_sink::write_atomically(
            &directory.join(format!("{name}.summary.tsv")),
            |mut file| Ok(file.write_all(summary.as_bytes())?),
        )?;
        triples(&directory.join(format!("{name}.links.zst")), &self.links)?;
        triples(&directory.join(format!("{name}.ties.zst")), &self.ties)
    }
}

fn ordered(first: u32, second: u32) -> (u32, u32) {
    (first.min(second), first.max(second))
}

fn triples(path: &Path, pairs: &HashMap<(u32, u32), u32>) -> Result<()> {
    let mut rows = pairs.iter().collect::<Vec<_>>();
    rows.sort_unstable();
    let bytes = rows
        .into_iter()
        .flat_map(|(&(first, second), &count)| [first, second, count])
        .flat_map(u32::to_le_bytes)
        .collect::<Vec<_>>();
    packed(path, |sink| Ok(sink.write_all(&bytes)?))
}
