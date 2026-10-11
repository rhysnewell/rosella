use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use log::info;
use ndarray::Array2;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::{
    self as sam, alignment::record::Cigar as _, header::record::value::map::header::tag::SORT_ORDER,
};
use rayon::prelude::*;

use super::coverage_table::CoverageTable;
use crate::learned::windows::{WINDOW, f16_bits};
use crate::progress::{self, Stage};

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];
pub const CONTIGS: &str = "contigs.tsv";

pub struct Sample {
    pub name: String,
    pub path: PathBuf,
}

impl Sample {
    pub fn windows(directory: &Path, name: &str) -> PathBuf {
        directory.join(format!("{name}.win.f16"))
    }

    fn table(&self, directory: &Path) -> PathBuf {
        directory.join(format!("{}.depth.tsv", self.name))
    }

    fn kept(&self, directory: &Path) -> Result<bool> {
        let read = std::fs::metadata(&self.path)?.modified()?;
        for path in [Self::windows(directory, &self.name), self.table(directory)] {
            match std::fs::metadata(path) {
                Ok(kept) if kept.modified()? >= read => {}
                _ => return Ok(false),
            }
        }
        Ok(true)
    }
}

#[derive(PartialEq)]
struct Contigs {
    names: Vec<String>,
    lengths: Vec<usize>,
}

// A BGZF file ends in a fixed empty block, so a BAM cut short by a killed sort is caught before
// its records are read as a complete sample.
pub fn complete(path: &Path) -> Result<bool> {
    let mut file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    if file.metadata()?.len() < BGZF_EOF.len() as u64 {
        return Ok(false);
    }
    file.seek(SeekFrom::End(-(BGZF_EOF.len() as i64)))?;
    let mut tail = [0; BGZF_EOF.len()];
    file.read_exact(&mut tail)?;
    Ok(tail == BGZF_EOF)
}

// Every primary mapped alignment counts over its whole reference span, with no identity filter
// and no end trimming, as the measured window depth was built.
pub fn read(samples: &[Sample], directory: &Path, threads: usize) -> Result<CoverageTable> {
    std::fs::create_dir_all(directory)?;
    let contigs = Contigs::of(&header_of(&samples[0].path)?);
    contigs.keep(directory)?;
    let mut fresh = Vec::new();
    for sample in samples {
        if !sample.kept(directory)? {
            if !complete(&sample.path)? {
                bail!(
                    "{} lacks the BGZF end block, so it was cut short",
                    sample.path.display()
                );
            }
            fresh.push(sample);
        }
    }
    if fresh.len() < samples.len() {
        info!(
            "Reusing the depth of {} samples kept in {}.",
            samples.len() - fresh.len(),
            directory.display()
        );
    }

    let readers = fresh.len().clamp(1, threads.max(1));
    let workers = threads / readers;
    let bar = progress::stepped(Stage::ReadingDepth, fresh.len() as u64);
    bar.set_message(format!("{} BAMs", fresh.len()));
    rayon::ThreadPoolBuilder::new()
        .num_threads(readers)
        .build()?
        .install(|| {
            fresh.par_iter().try_for_each(|sample| {
                let done = depth_of(sample, directory, &contigs, workers)
                    .with_context(|| format!("reading {}", sample.path.display()));
                bar.inc(1);
                done
            })
        })?;
    progress::done(&bar);

    let mut tables = samples
        .iter()
        .map(|sample| CoverageTable::from_any_file(sample.table(directory)));
    let mut table = tables.next().context("no BAMs to read")??;
    for next in tables {
        table.merge(next?)?;
    }
    Ok(table)
}

// Samples are taken in name order, so a directory filled one sample at a time reads the same
// however it was filled.
pub fn kept(directory: &Path) -> Result<CoverageTable> {
    let mut tables = std::fs::read_dir(directory)
        .with_context(|| format!("reading {}", directory.display()))?
        .map(|entry| Ok(entry?.path()))
        .collect::<Result<Vec<_>>>()?;
    // macOS leaves a hidden ._ twin beside each file on drives without extended attributes.
    tables.retain(|path| {
        path.file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.ends_with(".depth.tsv") && !name.starts_with('.'))
    });
    tables.sort();
    let mut tables = tables.iter().map(CoverageTable::from_any_file);
    let mut table = tables
        .next()
        .with_context(|| format!("{} holds no sample depth", directory.display()))??;
    for next in tables {
        table.merge(next?)?;
    }
    Ok(table)
}

fn header_of(path: &Path) -> Result<sam::Header> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    bam::io::Reader::new(BufReader::new(file))
        .read_header()
        .with_context(|| format!("reading the header of {}", path.display()))
}

impl Contigs {
    fn of(header: &sam::Header) -> Self {
        let (names, lengths) = header
            .reference_sequences()
            .iter()
            .map(|(name, sequence)| (name.to_string(), sequence.length().get()))
            .unzip();
        Self { names, lengths }
    }

    fn keep(&self, directory: &Path) -> Result<()> {
        let path = directory.join(CONTIGS);
        if path.is_file() {
            let (mut names, mut lengths) = (Vec::new(), Vec::new());
            for line in BufReader::new(File::open(&path)?).lines().skip(1) {
                let line = line?;
                let (name, length) = line.split_once('\t').context("a contig has no length")?;
                names.push(name.to_string());
                lengths.push(length.parse()?);
            }
            if *self != (Self { names, lengths }) {
                bail!(
                    "{} holds the depth of other contigs than these BAMs were mapped to, so \
                     remove it or name another --depth-directory",
                    directory.display()
                );
            }
            return Ok(());
        }
        let rows = self
            .names
            .iter()
            .zip(&self.lengths)
            .map(|(name, length)| format!("{name}\t{length}\n"))
            .collect::<String>();
        crate::report_sink::write_atomically(&path, |mut file| {
            Ok(file.write_all(format!("contig\tlength\n{rows}").as_bytes())?)
        })
    }
}

fn depth_of(sample: &Sample, directory: &Path, contigs: &Contigs, workers: usize) -> Result<()> {
    let file = BufReader::new(File::open(&sample.path)?);
    let windows = Sample::windows(directory, &sample.name);
    let partial = windows.with_extension("f16.partial");
    let sink = BufWriter::new(File::create(&partial)?);
    let moments = match NonZero::new(workers).filter(|workers| workers.get() > 1) {
        Some(workers) => tally(
            bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(
                workers, file,
            )),
            sample,
            contigs,
            sink,
        ),
        None => tally(bam::io::Reader::new(file), sample, contigs, sink),
    }?;
    std::fs::rename(&partial, &windows)?;

    let mut table = Array2::zeros((contigs.names.len(), 2));
    for (row, (mean, variance)) in moments.into_iter().enumerate() {
        table[[row, 0]] = mean;
        table[[row, 1]] = variance;
    }
    CoverageTable::new(
        contigs.names.clone(),
        contigs.lengths.clone(),
        vec![sample.name.clone()],
        table,
    )
    .write(sample.table(directory))
}

fn tally<R: Read, W: Write>(
    mut reader: bam::io::Reader<R>,
    sample: &Sample,
    contigs: &Contigs,
    sink: W,
) -> Result<Vec<(f64, f64)>> {
    let path = sample.path.display();
    let header = reader.read_header()?;
    if let Some(order) = header
        .header()
        .and_then(|map| map.other_fields().get(&SORT_ORDER))
        .filter(|order| order.as_slice() != b"coordinate")
    {
        bail!("{path} is sorted by {order}, sort it by coordinate first");
    }
    if Contigs::of(&header) != *contigs {
        bail!("{path} was mapped to other contigs than the first BAM");
    }
    let mut pileup = Pileup::new(&contigs.lengths, sink);
    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        let (Some(contig), Some(start)) = (
            record.reference_sequence_id().transpose()?,
            record.alignment_start().transpose()?,
        ) else {
            continue;
        };
        let span = record.cigar().alignment_span()?;
        if !pileup.add(contig, usize::from(start) - 1, span)? {
            bail!("{path} is not sorted by coordinate, sort it first");
        }
    }
    pileup.finish()
}

struct Pileup<'a, W: Write> {
    lengths: &'a [usize],
    sink: W,
    change: Vec<i32>,
    current: Option<usize>,
    written: usize,
    moments: Vec<(f64, f64)>,
}

impl<'a, W: Write> Pileup<'a, W> {
    fn new(lengths: &'a [usize], sink: W) -> Self {
        Self {
            lengths,
            sink,
            change: vec![0; lengths.iter().max().map_or(0, |longest| longest + 1)],
            current: None,
            written: 0,
            moments: vec![(0.0, 0.0); lengths.len()],
        }
    }

    fn add(&mut self, contig: usize, start: usize, span: usize) -> Result<bool> {
        if self.current != Some(contig) {
            if self.current.is_some_and(|current| contig < current) {
                return Ok(false);
            }
            self.close()?;
            self.current = Some(contig);
        }
        let end = (start + span).min(self.lengths[contig]);
        if start < end {
            self.change[start] += 1;
            self.change[end] -= 1;
        }
        Ok(true)
    }

    fn finish(mut self) -> Result<Vec<(f64, f64)>> {
        self.close()?;
        self.blank_to(self.lengths.len())?;
        self.sink.flush()?;
        Ok(self.moments)
    }

    fn blank_to(&mut self, contig: usize) -> Result<()> {
        let windows = self.lengths[self.written..contig]
            .iter()
            .map(|length| length.div_ceil(WINDOW))
            .sum::<usize>();
        self.sink.write_all(&vec![0; 2 * windows])?;
        self.written = contig;
        Ok(())
    }

    fn close(&mut self) -> Result<()> {
        let Some(contig) = self.current.take() else {
            return Ok(());
        };
        self.blank_to(contig)?;
        let length = self.lengths[contig];
        let (mut depth, mut sum, mut squares, mut window_sum) = (0i64, 0.0, 0.0, 0i64);
        let mut windows = Vec::with_capacity(2 * length.div_ceil(WINDOW));
        for (position, delta) in self.change[..length].iter_mut().enumerate() {
            depth += i64::from(*delta);
            *delta = 0;
            let value = depth as f64;
            sum += value;
            squares += value * value;
            window_sum += depth;
            if (position + 1) % WINDOW == 0 || position + 1 == length {
                let width = position % WINDOW + 1;
                windows.extend(f16_bits(window_sum as f32 / width as f32).to_le_bytes());
                window_sum = 0;
            }
        }
        self.change[length] = 0;
        self.sink.write_all(&windows)?;
        self.written = contig + 1;
        let mean = sum / length as f64;
        self.moments[contig] = (mean, (squares / length as f64 - mean * mean).max(0.0));
        Ok(())
    }
}
