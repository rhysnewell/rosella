use std::collections::HashMap;
use std::fs::File;
use std::hash::Hasher;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use log::info;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam as sam;

use super::coverage_table::CoverageTable;
use super::depth_file::{self, Contigs, Depth};
use super::tie_split;
use crate::digest::Fnv1a;
use crate::learned::windows::{WINDOW, f16_bits};
use crate::progress::{self, Stage};

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];
const EDGE: u64 = 1 << 16;

pub struct Sample {
    pub name: String,
    pub path: PathBuf,
}

impl Sample {
    // rsync -a carries a BAM's mtime with it, so a BAM copied over the one that made the depth is
    // told apart by its bytes.
    fn kept(&self, directory: &Path, key: &str) -> bool {
        depth_file::path(directory, &self.name).is_file()
            && std::fs::read_to_string(key_path(directory, &self.name))
                .is_ok_and(|kept| kept == key)
    }
}

fn key_path(directory: &Path, name: &str) -> PathBuf {
    directory.join(format!("{name}.bam.key"))
}

// The head holds the mapper's settings and the tail the last hits, so with the length they tell
// one mapping from another without reading the whole BAM.
fn key(path: &Path) -> Result<String> {
    let mut file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let length = file.metadata()?.len();
    let mut edges = Vec::with_capacity(2 * EDGE as usize);
    (&mut file).take(EDGE).read_to_end(&mut edges)?;
    file.seek(SeekFrom::Start(length.saturating_sub(EDGE)))?;
    file.take(EDGE).read_to_end(&mut edges)?;
    let mut hash = Fnv1a::default();
    hash.write(&length.to_le_bytes());
    hash.write(&edges);
    Ok(format!("{:016x}", hash.finish()))
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

// No identity filter and no end trimming, as the measured window depth was built. Every BAM
// holds its hits in memory for the tie split, so they are read one at a time.
pub fn read(samples: &[Sample], directory: &Path, threads: usize) -> Result<CoverageTable> {
    std::fs::create_dir_all(directory)?;
    let contigs = contigs_of(&header_of(&samples[0].path)?);
    contigs.keep(directory)?;
    let index = contigs
        .names
        .iter()
        .enumerate()
        .map(|(at, name)| (name.as_str(), at as u32))
        .collect::<HashMap<_, _>>();
    let mut fresh = Vec::new();
    for sample in samples {
        let key = key(&sample.path)?;
        if !sample.kept(directory, &key) {
            if !complete(&sample.path)? {
                bail!(
                    "{} is cut short and lacks the BGZF end block",
                    sample.path.display()
                );
            }
            fresh.push((sample, key));
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
    let reading = Reading {
        directory,
        contigs: &contigs,
        index: &index,
    };
    for (sample, key) in fresh {
        reading.depth_of(sample, &key, in_read_order(&sample.path)?, threads)?;
        bar.inc(1);
    }
    progress::done(&bar);

    let names = samples
        .iter()
        .map(|sample| sample.name.clone())
        .collect::<Vec<_>>();
    depth_file::merged(directory, &contigs, &names)
}

pub fn header_of(path: &Path) -> Result<sam::Header> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    bam::io::Reader::new(BufReader::new(file))
        .read_header()
        .with_context(|| format!("reading the header of {}", path.display()))
}

// The order sits on the header's first line, so it is read without parsing every contig name.
fn in_read_order(path: &Path) -> Result<bool> {
    let opened = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut reader = bgzf::io::Reader::new(opened);
    let mut head = [0; 8];
    reader
        .read_exact(&mut head)
        .with_context(|| format!("reading the header of {}", path.display()))?;
    if &head[..4] != b"BAM\x01" {
        bail!("{} is not a BAM", path.display());
    }
    let text = u32::from_le_bytes([head[4], head[5], head[6], head[7]]);
    let mut line = Vec::new();
    reader.take(u64::from(text)).read_until(b'\n', &mut line)?;
    Ok(line.starts_with(b"@HD")
        && line
            .split(|byte| matches!(byte, b'\t' | b'\n' | b'\r'))
            .any(|field| field == b"SO:queryname" || field == b"GO:query"))
}

fn contigs_of(header: &sam::Header) -> Contigs {
    let (names, lengths) = header
        .reference_sequences()
        .iter()
        .map(|(name, sequence)| (name.to_string(), sequence.length().get()))
        .unzip();
    Contigs { names, lengths }
}

struct Reading<'a> {
    directory: &'a Path,
    contigs: &'a Contigs,
    index: &'a HashMap<&'a str, u32>,
}

impl Reading<'_> {
    fn depth_of(&self, sample: &Sample, key: &str, grouped: bool, workers: usize) -> Result<()> {
        let open = || self.opened(&sample.path, workers);
        let (pileup, stats) = match grouped {
            true => open().and_then(|reader| tie_split::tally(reader, self.contigs, self.index)),
            false => tie_split::sorted(open, self.contigs, self.index),
        }
        .with_context(|| format!("reading {}", sample.path.display()))?;
        stats.write(self.directory, &sample.name)?;
        let (moments, windows) = pileup.finish();
        let depth = Depth {
            means: moments.iter().map(|moment| moment.0).collect(),
            variances: moments.iter().map(|moment| moment.1).collect(),
            reads: stats.reads,
        };
        depth_file::write(
            &depth_file::path(self.directory, &sample.name),
            &depth,
            &windows,
        )?;
        crate::report_sink::write_atomically(&key_path(self.directory, &sample.name), |mut file| {
            Ok(file.write_all(key.as_bytes())?)
        })
    }

    fn opened(&self, path: &Path, workers: usize) -> Result<bam::io::Reader<Box<dyn Read>>> {
        let file = BufReader::new(File::open(path)?);
        let blocks: Box<dyn Read> = match NonZero::new(workers).filter(|workers| workers.get() > 1)
        {
            Some(workers) => Box::new(bgzf::io::MultithreadedReader::with_worker_count(
                workers, file,
            )),
            None => Box::new(bgzf::io::Reader::new(file)),
        };
        let mut reader = bam::io::Reader::from(blocks);
        let header = reader.read_header()?;
        let same = header
            .reference_sequences()
            .iter()
            .map(|(name, sequence)| (name.as_slice(), sequence.length().get()))
            .eq(self
                .contigs
                .names
                .iter()
                .map(String::as_bytes)
                .zip(self.contigs.lengths.iter().copied()));
        if !same {
            bail!(
                "{} was mapped to other contigs than the first BAM",
                path.display()
            );
        }
        Ok(reader)
    }
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

    // Contigs come in order, since the tie split lays every hit out contig by contig.
    pub(super) fn add(&mut self, contig: usize, start: usize, span: usize, weight: f64) {
        if self.current != Some(contig) {
            self.close();
            self.current = Some(contig);
        }
        let end = (start + span).min(self.lengths[contig]);
        if start < end {
            self.change[start] += weight;
            self.change[end] -= weight;
        }
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
