use std::ops::Range;
use std::sync::mpsc::sync_channel;
use std::thread;
use std::time::{Duration, Instant};

use anyhow::Result;
use frugal::api::{MetaPredictor, ProdigalConfig, Strand};
use log::debug;
use needletail::parse_fastx_file;
use rayon::prelude::*;

use crate::pool;

pub struct Orf {
    pub contig: usize,
    pub bases: usize,
    pub partial: bool,
    pub protein: String,
}

/// Contigs are called two chunks at a time so the assembly is never held whole beside the
/// proteins it produces, which was the run's memory peak on a multi-sample assembly.
const CHUNK_BASES: usize = 64 << 20;

#[derive(Default)]
struct Chunk {
    names: Vec<String>,
    lengths: Vec<usize>,
    slots: Vec<usize>,
    held: Vec<Vec<u8>>,
}

#[derive(Debug, Default)]
pub struct Walked {
    pub names: Vec<String>,
    pub lengths: Vec<usize>,
}

type Chunks = std::sync::mpsc::SyncSender<Result<Chunk>>;

// Contigs at or above `band.end` are named but not called, so an annotation held above a floor
// extends below it without calling a contig twice.
pub fn call_over<F>(assembly: &str, band: Range<usize>, mut batch: F) -> Result<Walked>
where
    F: FnMut(Vec<Orf>) -> Result<()>,
{
    let predictor =
        MetaPredictor::with_config_and_thread_pool(ProdigalConfig::default(), pool::get())
            .map_err(|error| anyhow!("gene finder: {error:?}"))?;

    let (sender, receiver) = sync_channel::<Result<Chunk>>(1);
    let held_assembly = assembly.to_string();
    let reader = thread::spawn(move || {
        if let Err(error) = read_chunks(&held_assembly, band, &sender) {
            let _ = sender.send(Err(error));
        }
    });

    let progress = crate::progress::spinning(crate::progress::Stage::CallingGenes);
    let mut walked = Walked::default();
    let mut outcome = Ok(());
    let mut spent = Spent::default();
    loop {
        let blocked = Instant::now();
        let Ok(chunk) = receiver.recv() else {
            break;
        };
        spent.waiting += blocked.elapsed();
        let mut chunk = match chunk {
            Ok(chunk) => chunk,
            Err(error) => {
                outcome = Err(error);
                break;
            }
        };
        walked.names.append(&mut chunk.names);
        walked.lengths.append(&mut chunk.lengths);
        progress.set_message(format!("{} contigs", walked.names.len()));
        let started = Instant::now();
        let called = call(&predictor, &chunk.slots, &chunk.held);
        spent.calling += started.elapsed();
        let started = Instant::now();
        outcome = called.and_then(&mut batch);
        spent.writing += started.elapsed();
        if outcome.is_err() {
            break;
        }
    }
    progress.finish_and_clear();
    drop(receiver);
    reader
        .join()
        .map_err(|_| anyhow!("the assembly reader panicked"))?;
    spent.report();
    outcome.map(|()| walked)
}

#[derive(Default)]
struct Spent {
    waiting: Duration,
    calling: Duration,
    writing: Duration,
}

impl Spent {
    fn report(&self) {
        debug!(
            "gene calling spent {:.1}s waiting on the reader, {:.1}s calling, {:.1}s in the callback",
            self.waiting.as_secs_f64(),
            self.calling.as_secs_f64(),
            self.writing.as_secs_f64()
        );
    }
}

fn read_chunks(assembly: &str, band: Range<usize>, sender: &Chunks) -> Result<()> {
    let mut reader = parse_fastx_file(assembly)?;
    let mut chunk = Chunk::default();
    let mut seen = 0usize;
    let mut bases = 0usize;
    while let Some(record) = reader.next() {
        let record = record?;
        let sequence = record.seq();
        if sequence.len() < band.start {
            continue;
        }
        chunk.names.push(crate::contig_id(record.id())?.to_string());
        chunk.lengths.push(sequence.len());
        seen += 1;
        if sequence.len() >= band.end {
            continue;
        }
        bases += sequence.len();
        chunk.held.push(sequence.to_vec());
        chunk.slots.push(seen - 1);
        if bases >= CHUNK_BASES {
            if sender.send(Ok(std::mem::take(&mut chunk))).is_err() {
                return Ok(());
            }
            bases = 0;
        }
    }
    if !chunk.names.is_empty() {
        let _ = sender.send(Ok(chunk));
    }
    Ok(())
}

/// Rayon splits the batch into contiguous ranges and folds each one to completion, so any
/// length-ordered array hands a single worker every long contig while the rest sleep.
const SPREAD_SEED: u64 = 0x2545_f491_4f6c_dd1d;

fn spread(count: usize) -> Vec<usize> {
    let mut order = (0..count).collect::<Vec<_>>();
    let mut state = SPREAD_SEED;
    for at in (1..count).rev() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        order.swap(at, (state % (at as u64 + 1)) as usize);
    }
    order
}

fn call(predictor: &MetaPredictor, slots: &[usize], contigs: &[Vec<u8>]) -> Result<Vec<Orf>> {
    let order = spread(contigs.len());
    let sequences = order
        .iter()
        .map(|index| contigs[*index].as_slice())
        .collect::<Vec<_>>();
    let batches = predictor
        .predict_batch(&sequences)
        .map_err(|error| anyhow!("gene finder: {error:?}"))?;

    let mut called = vec![Vec::new(); contigs.len()];
    for (slot, genes) in order.iter().zip(batches) {
        called[*slot] = genes;
    }

    let translated = contigs
        .par_iter()
        .zip(called.into_par_iter())
        .enumerate()
        .map(|(at, (sequence, genes))| {
            genes
                .into_iter()
                .map(|gene| {
                    let coding = &sequence[gene.begin - 1..gene.end];
                    let protein = match gene.strand {
                        Strand::Reverse => translate_reverse(coding, !gene.partial.1),
                        _ => translate(coding, !gene.partial.0),
                    };
                    Orf {
                        contig: slots[at],
                        bases: gene.end - gene.begin + 1,
                        partial: gene.partial.0 || gene.partial.1,
                        protein,
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    Ok(translated.into_iter().flatten().collect())
}

/// NCBI table 11 differs from the standard code only in which codons may initiate, so the
/// residue mapping is the standard one and the start is rewritten separately.
const CODONS: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

fn base(byte: u8) -> Option<usize> {
    match byte.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

fn complement(byte: u8) -> Option<usize> {
    base(byte).map(|index| (index + 2) % 4)
}

fn finish(mut protein: String, complete_start: bool) -> String {
    if protein.ends_with('*') {
        protein.pop();
    }
    if complete_start && !protein.is_empty() {
        protein.replace_range(0..1, "M");
    }
    protein
}

fn residue(first: Option<usize>, second: Option<usize>, third: Option<usize>) -> char {
    match (first, second, third) {
        (Some(a), Some(b), Some(c)) => CODONS[a * 16 + b * 4 + c] as char,
        _ => 'X',
    }
}

pub fn translate(coding: &[u8], complete_start: bool) -> String {
    let mut protein = String::with_capacity(coding.len() / 3);
    for codon in coding.chunks_exact(3) {
        protein.push(residue(base(codon[0]), base(codon[1]), base(codon[2])));
    }
    finish(protein, complete_start)
}

pub fn translate_reverse(coding: &[u8], complete_start: bool) -> String {
    let mut protein = String::with_capacity(coding.len() / 3);
    let mut end = coding.len();
    while end >= 3 {
        protein.push(residue(
            complement(coding[end - 1]),
            complement(coding[end - 2]),
            complement(coding[end - 3]),
        ));
        end -= 3;
    }
    finish(protein, complete_start)
}
