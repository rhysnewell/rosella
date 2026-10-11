use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use log::info;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::{
    self as sam, alignment::record::Cigar as _, header::record::value::map::header::tag::SORT_ORDER,
};
use rayon::prelude::*;

use super::coverage_table::CoverageTable;
use super::depth_file::{self, Contigs, Depth};
use crate::learned::windows::{WINDOW, f16_bits};
use crate::progress::{self, Stage};

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

pub struct Sample {
    pub name: String,
    pub path: PathBuf,
}

impl Sample {
    fn kept(&self, directory: &Path) -> Result<bool> {
        let read = std::fs::metadata(&self.path)?.modified()?;
        Ok(std::fs::metadata(depth_file::path(directory, &self.name))
            .and_then(|kept| kept.modified())
            .is_ok_and(|kept| kept >= read))
    }
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
    let contigs = contigs_of(&header_of(&samples[0].path)?);
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

    let names = samples
        .iter()
        .map(|sample| sample.name.clone())
        .collect::<Vec<_>>();
    depth_file::merged(directory, &contigs, &names)
}

fn header_of(path: &Path) -> Result<sam::Header> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    bam::io::Reader::new(BufReader::new(file))
        .read_header()
        .with_context(|| format!("reading the header of {}", path.display()))
}

fn contigs_of(header: &sam::Header) -> Contigs {
    let (names, lengths) = header
        .reference_sequences()
        .iter()
        .map(|(name, sequence)| (name.to_string(), sequence.length().get()))
        .unzip();
    Contigs { names, lengths }
}

fn depth_of(sample: &Sample, directory: &Path, contigs: &Contigs, workers: usize) -> Result<()> {
    let file = BufReader::new(File::open(&sample.path)?);
    let pileup = match NonZero::new(workers).filter(|workers| workers.get() > 1) {
        Some(workers) => tally(
            bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(
                workers, file,
            )),
            sample,
            contigs,
        ),
        None => tally(bam::io::Reader::new(file), sample, contigs),
    }?;
    let depth = Depth {
        means: pileup.moments.iter().map(|moment| moment.0).collect(),
        variances: pileup.moments.iter().map(|moment| moment.1).collect(),
        reads: pileup.reads,
    };
    depth_file::write(
        &depth_file::path(directory, &sample.name),
        &depth,
        &pileup.windows,
    )
}

fn tally<'a, R: Read>(
    mut reader: bam::io::Reader<R>,
    sample: &Sample,
    contigs: &'a Contigs,
) -> Result<Pileup<'a>> {
    let path = sample.path.display();
    let header = reader.read_header()?;
    if let Some(order) = header
        .header()
        .and_then(|map| map.other_fields().get(&SORT_ORDER))
        .filter(|order| order.as_slice() != b"coordinate")
    {
        bail!("{path} is sorted by {order}, sort it by coordinate first");
    }
    if contigs_of(&header) != *contigs {
        bail!("{path} was mapped to other contigs than the first BAM");
    }
    let mut pileup = Pileup::new(contigs);
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
        if !pileup.add(contig, usize::from(start) - 1, span) {
            bail!("{path} is not sorted by coordinate, sort it first");
        }
    }
    pileup.close();
    pileup.blank_to(contigs.lengths.len());
    Ok(pileup)
}

// The windows are held in memory, two bytes per 100 bp, so the per-contig columns can be
// written ahead of them.
struct Pileup<'a> {
    lengths: &'a [usize],
    change: Vec<i32>,
    current: Option<usize>,
    written: usize,
    windows: Vec<u8>,
    moments: Vec<(f64, f64)>,
    reads: Vec<u32>,
}

impl<'a> Pileup<'a> {
    fn new(contigs: &'a Contigs) -> Self {
        let lengths = &contigs.lengths;
        Self {
            lengths,
            change: vec![0; lengths.iter().max().map_or(0, |longest| longest + 1)],
            current: None,
            written: 0,
            windows: Vec::with_capacity(2 * contigs.windows()),
            moments: vec![(0.0, 0.0); lengths.len()],
            reads: vec![0; lengths.len()],
        }
    }

    fn add(&mut self, contig: usize, start: usize, span: usize) -> bool {
        if self.current != Some(contig) {
            if self.current.is_some_and(|current| contig < current) {
                return false;
            }
            self.close();
            self.current = Some(contig);
        }
        self.reads[contig] += 1;
        let end = (start + span).min(self.lengths[contig]);
        if start < end {
            self.change[start] += 1;
            self.change[end] -= 1;
        }
        true
    }

    fn blank_to(&mut self, contig: usize) {
        let windows = self.lengths[self.written..contig]
            .iter()
            .map(|length| length.div_ceil(WINDOW))
            .sum::<usize>();
        self.windows.resize(self.windows.len() + 2 * windows, 0);
        self.written = contig;
    }

    fn close(&mut self) {
        let Some(contig) = self.current.take() else {
            return;
        };
        self.blank_to(contig);
        let length = self.lengths[contig];
        let (mut depth, mut sum, mut squares, mut window_sum) = (0i64, 0.0, 0.0, 0i64);
        for (position, delta) in self.change[..length].iter_mut().enumerate() {
            depth += i64::from(*delta);
            *delta = 0;
            let value = depth as f64;
            sum += value;
            squares += value * value;
            window_sum += depth;
            if (position + 1) % WINDOW == 0 || position + 1 == length {
                let width = position % WINDOW + 1;
                self.windows
                    .extend(f16_bits(window_sum as f32 / width as f32).to_le_bytes());
                window_sum = 0;
            }
        }
        self.change[length] = 0;
        self.written = contig + 1;
        let mean = sum / length as f64;
        self.moments[contig] = (mean, (squares / length as f64 - mean * mean).max(0.0));
    }
}
