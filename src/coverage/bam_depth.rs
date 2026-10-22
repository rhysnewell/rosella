use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use log::info;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::{
    self as sam,
    alignment::record::Cigar as _,
    header::record::value::map::header::tag::{GROUP_ORDER, SORT_ORDER},
};
use rayon::prelude::*;

use super::coverage_table::CoverageTable;
use super::depth_file::{self, Contigs, Depth};
use super::read_stats::ReadStats;
use super::tie_split;
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

// No identity filter and no end trimming, as the measured window depth was built. A BAM in
// mapper order holds its tied reads in memory, so those are read one at a time.
pub fn read(samples: &[Sample], directory: &Path, threads: usize) -> Result<CoverageTable> {
    std::fs::create_dir_all(directory)?;
    let contigs = contigs_of(&header_of(&samples[0].path)?);
    contigs.keep(directory)?;
    let mut fresh = Vec::new();
    for sample in samples {
        if !sample.kept(directory)? {
            if !complete(&sample.path)? {
                bail!(
                    "{} is cut short and lacks the BGZF end block",
                    sample.path.display()
                );
            }
            fresh.push(sample);
        }
    }
    if fresh.len() < samples.len() {
        info!(
            "Depth reused: {} samples from {}",
            samples.len() - fresh.len(),
            directory.display()
        );
    }

    let bar = progress::stepped(Stage::ReadingDepth, fresh.len() as u64);
    bar.set_message(format!("{} BAMs", fresh.len()));
    let (mut sorted, mut grouped) = (Vec::new(), Vec::new());
    for sample in fresh {
        match in_read_order(&header_of(&sample.path)?) {
            true => grouped.push(sample),
            false => sorted.push(sample),
        }
    }
    let readers = sorted.len().clamp(1, threads.max(1));
    rayon::ThreadPoolBuilder::new()
        .num_threads(readers)
        .build()?
        .install(|| {
            sorted.par_iter().try_for_each(|sample| {
                let done = depth_of(sample, directory, &contigs, threads / readers);
                bar.inc(1);
                done
            })
        })?;
    for sample in grouped {
        depth_of(sample, directory, &contigs, threads)?;
        bar.inc(1);
    }
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

fn in_read_order(header: &sam::Header) -> bool {
    let field = |tag| {
        header
            .header()
            .and_then(|map| map.other_fields().get(&tag))
            .map(|value| value.as_slice())
    };
    field(SORT_ORDER) == Some(b"queryname") || field(GROUP_ORDER) == Some(b"query")
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
    let (pileup, stats) = match NonZero::new(workers).filter(|workers| workers.get() > 1) {
        Some(workers) => tally(
            bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(
                workers, file,
            )),
            sample,
            contigs,
        ),
        None => tally(bam::io::Reader::new(file), sample, contigs),
    }
    .with_context(|| format!("reading {}", sample.path.display()))?;
    stats.write(directory, &sample.name)?;
    let (moments, windows) = pileup.finish();
    let depth = Depth {
        means: moments.iter().map(|moment| moment.0).collect(),
        variances: moments.iter().map(|moment| moment.1).collect(),
        reads: stats.reads,
    };
    depth_file::write(&depth_file::path(directory, &sample.name), &depth, &windows)
}

fn tally<'a, R: Read>(
    mut reader: bam::io::Reader<R>,
    sample: &Sample,
    contigs: &'a Contigs,
) -> Result<(Pileup<'a>, ReadStats<'a>)> {
    let path = sample.path.display();
    let header = reader.read_header()?;
    if contigs_of(&header) != *contigs {
        bail!("{path} was mapped to other contigs than the first BAM");
    }
    if in_read_order(&header) {
        return tie_split::tally(reader, contigs);
    }
    if let Some(order) = header
        .header()
        .and_then(|map| map.other_fields().get(&SORT_ORDER))
        .filter(|order| order.as_slice() != b"coordinate")
    {
        bail!("{path} is sorted by {order}. Sort it by coordinate first");
    }
    let mut pileup = Pileup::new(contigs);
    let mut stats = ReadStats::new(contigs);
    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        stats.note(&record)?;
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
        if !pileup.add(contig, usize::from(start) - 1, span, 1.0) {
            bail!("{path} is not sorted by coordinate. Sort it first");
        }
    }
    Ok((pileup, stats))
}

// The windows are held in memory, two bytes per 100 bp, so the per-contig columns can be
// written ahead of them.
pub(super) struct Pileup<'a> {
    lengths: &'a [usize],
    change: Vec<f64>,
    current: Option<usize>,
    written: usize,
    windows: Vec<u8>,
    moments: Vec<(f64, f64)>,
}

impl<'a> Pileup<'a> {
    pub(super) fn new(contigs: &'a Contigs) -> Self {
        let lengths = &contigs.lengths;
        Self {
            lengths,
            change: vec![0.0; lengths.iter().max().map_or(0, |longest| longest + 1)],
            current: None,
            written: 0,
            windows: Vec::with_capacity(2 * contigs.windows()),
            moments: vec![(0.0, 0.0); lengths.len()],
        }
    }

    pub(super) fn add(&mut self, contig: usize, start: usize, span: usize, weight: f64) -> bool {
        if self.current != Some(contig) {
            if self.current.is_some_and(|current| contig < current) {
                return false;
            }
            self.close();
            self.current = Some(contig);
        }
        let end = (start + span).min(self.lengths[contig]);
        if start < end {
            self.change[start] += weight;
            self.change[end] -= weight;
        }
        true
    }

    fn finish(mut self) -> (Vec<(f64, f64)>, Vec<u8>) {
        self.close();
        self.blank_to(self.lengths.len());
        (self.moments, self.windows)
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
        let (mut depth, mut sum, mut squares, mut window_sum) = (0.0, 0.0, 0.0, 0.0);
        for (position, delta) in self.change[..length].iter_mut().enumerate() {
            depth += *delta;
            *delta = 0.0;
            // Weighted reads leave a rounding residue where they end, which must not read as depth.
            let value = f64::max(depth, 0.0);
            sum += value;
            squares += value * value;
            window_sum += value;
            if (position + 1) % WINDOW == 0 || position + 1 == length {
                let width = position % WINDOW + 1;
                self.windows
                    .extend(f16_bits(window_sum as f32 / width as f32).to_le_bytes());
                window_sum = 0.0;
            }
        }
        self.change[length] = 0.0;
        self.written = contig + 1;
        let mean = sum / length as f64;
        self.moments[contig] = (mean, (squares / length as f64 - mean * mean).max(0.0));
    }
}
