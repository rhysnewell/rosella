use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

use anyhow::{Context, Result};
use log::info;
use noodles_bam as bam;
use noodles_sam::alignment::record::{
    Cigar as _,
    data::field::{Tag, Value},
};
use rayon::prelude::*;

use super::bridges::{Piece, bridges};
use super::depth_file::{self, Contigs, packed};

pub(super) struct ReadStats<'a> {
    names: &'a HashMap<&'a str, usize>,
    lengths: &'a [usize],
    pub(super) reads: Vec<u32>,
    seen: u64,
    mapped: u64,
    tied: u64,
    read_bases: u64,
    mapped_bases: u64,
    long: bool,
    bridges: HashMap<(u32, u32), u32>,
    ties: HashMap<(u32, u32), u32>,
}

impl<'a> ReadStats<'a> {
    pub(super) fn new(
        contigs: &'a Contigs,
        names: &'a HashMap<&'a str, usize>,
        long: bool,
    ) -> Self {
        Self {
            names,
            long,
            lengths: &contigs.lengths,
            reads: vec![0; contigs.names.len()],
            seen: 0,
            mapped: 0,
            tied: 0,
            read_bases: 0,
            mapped_bases: 0,
            bridges: HashMap::new(),
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
        if self.long
            && let Some(Value::String(others)) =
                record.data().get(&Tag::OTHER_ALIGNMENTS).transpose()?
            && let Some(start) = record.alignment_start().transpose()?
        {
            let ops = cigar
                .iter()
                .map(|op| op.map(|op| (op.kind(), op.len())))
                .collect::<std::io::Result<Vec<_>>>()?;
            let own = Piece::of(usize::from(start) - 1, ops, flags.is_reverse_complemented());
            for other in others.split(|byte| *byte == b';') {
                let Some((linked, piece)) = self.alignment(other) else {
                    continue;
                };
                if linked != contig
                    && bridges(&own, self.lengths[contig], &piece, self.lengths[linked])
                {
                    *self
                        .bridges
                        .entry(ordered(contig as u32, linked as u32))
                        .or_default() += 1;
                }
            }
        }
        Ok(())
    }

    fn alignment(&self, entry: &[u8]) -> Option<(usize, Piece)> {
        let mut fields = std::str::from_utf8(entry).ok()?.split(',');
        let linked = *self.names.get(fields.next()?)?;
        let start = fields.next()?.parse::<usize>().ok()?.checked_sub(1)?;
        let reverse = fields.next()? == "-";
        Some((linked, Piece::parsed(start, fields.next()?, reverse)?))
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
        let bridged = bridges_path(directory, name);
        match self.long {
            true => triples(&bridged, &self.bridges)?,
            // A short sample's bridges left by an older make would still be read as long reads.
            false => match std::fs::remove_file(&bridged) {
                Err(error) if error.kind() != std::io::ErrorKind::NotFound => Err(error)?,
                _ => {}
            },
        }
        triples(&directory.join(format!("{name}.ties.zst")), &self.ties)
    }
}

pub(super) fn bridges_path(directory: &Path, name: &str) -> std::path::PathBuf {
    directory.join(format!("{name}.bridges.zst"))
}

pub struct Placed {
    pub at: Vec<Option<usize>>,
    pub reads: Vec<u64>,
}

impl Placed {
    pub fn read(directory: &Path, samples: &[String], names: &[String]) -> Result<Self> {
        let contigs = Contigs::read(directory)?;
        let index = crate::rows::positions(contigs.names.iter().map(String::as_str));
        let at = names
            .iter()
            .map(|name| index.get(name.as_str()).copied())
            .collect();
        let mut reads = vec![0u64; contigs.names.len()];
        for sample in samples {
            let depth = depth_file::read(&depth_file::path(directory, sample), reads.len())?;
            reads
                .iter_mut()
                .zip(depth.reads)
                .for_each(|(sum, own)| *sum += u64::from(own));
        }
        Ok(Self { at, reads })
    }

    pub fn wanted(&self, rows: impl IntoIterator<Item = usize>) -> Vec<bool> {
        let mut wanted = vec![false; self.reads.len()];
        for row in rows {
            if let Some(contig) = self.at[row] {
                wanted[contig] = true;
            }
        }
        wanted
    }
}

fn pair_rows(directory: &Path, sample: &str, kind: &str) -> Result<Vec<[u32; 3]>> {
    let path = directory.join(format!("{sample}.{kind}.zst"));
    if !path.is_file() {
        return Ok(Vec::new());
    }
    let mut bytes = Vec::new();
    zstd::Decoder::new(File::open(&path)?)?
        .read_to_end(&mut bytes)
        .with_context(|| format!("reading {}", path.display()))?;
    let words = bytes
        .as_chunks::<4>()
        .0
        .iter()
        .map(|word| u32::from_le_bytes(*word))
        .collect::<Vec<_>>();
    Ok(words.as_chunks::<3>().0.to_vec())
}

// Reads each wanted contig ties with another contig, over every sample. A read tied among k
// contigs counts k - 1 times, which only runs high where many strains share sequence.
pub fn tied_reads(directory: &Path, samples: &[String], wanted: &[bool]) -> Result<Vec<u64>> {
    let per_sample = samples
        .par_iter()
        .map(|sample| {
            let mut tied = vec![0u64; wanted.len()];
            for [first, second, count] in pair_rows(directory, sample, "ties")? {
                for contig in [first, second] {
                    if wanted.get(contig as usize).copied().unwrap_or(false) {
                        tied[contig as usize] += u64::from(count);
                    }
                }
            }
            Ok(tied)
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(per_sample
        .into_iter()
        .reduce(|mut sum, tied| {
            sum.iter_mut()
                .zip(tied)
                .for_each(|(sum, tied)| *sum += tied);
            sum
        })
        .unwrap_or_else(|| vec![0; wanted.len()]))
}

// A depth directory from an older make has no summary, so the read lengths are all or nothing.
pub fn read_lengths(directory: &Path, samples: &[String]) -> Result<Option<Vec<f64>>> {
    let mut lengths = Vec::with_capacity(samples.len());
    for sample in samples {
        let path = directory.join(format!("{sample}.summary.tsv"));
        if !path.is_file() {
            return Ok(None);
        }
        let text = std::fs::read_to_string(&path)?;
        let mut lines = text.lines();
        let column = lines
            .next()
            .and_then(|header| {
                header
                    .split('\t')
                    .position(|field| field == "mean_read_length")
            })
            .with_context(|| format!("no mean_read_length in {}", path.display()))?;
        let length = lines
            .next()
            .and_then(|row| row.split('\t').nth(column))
            .and_then(|field| field.parse::<f64>().ok())
            .with_context(|| format!("no read length row in {}", path.display()))?;
        lengths.push(length);
    }
    Ok(Some(lengths))
}

pub fn bridged(directory: &Path, samples: &[String]) -> Result<HashMap<(u32, u32), u32>> {
    let mut pairs = HashMap::new();
    for sample in samples {
        for [first, second, count] in pair_rows(directory, sample, "bridges")? {
            *pairs.entry((first, second)).or_default() += count;
        }
    }
    Ok(pairs)
}

// A tied read's primary lands on one of its contigs, so a contig's own count holds only half of
// the reads it ties. Exact when each read ties between two contigs.
pub fn tie_share(tied: u64, reads: u64) -> f64 {
    match tied {
        0 => 0.0,
        _ => (tied as f64 / (reads as f64 + tied as f64 / 2.0)).min(1.0),
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
