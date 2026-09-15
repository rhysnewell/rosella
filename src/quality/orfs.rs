use std::sync::mpsc::sync_channel;
use std::thread;

use anyhow::Result;
use needletail::parse_fastx_file;
use frugal::api::{MetaPredictor, ProdigalConfig, Strand};
use rayon::prelude::*;

use crate::pool;

pub struct Orf {
    pub contig: usize,
    pub partial: bool,
    pub protein: String,
}

/// Contigs are called two chunks at a time so the assembly is never held whole beside the
/// proteins it produces, which was the run's memory peak on a multi-sample assembly.
const CHUNK_BASES: usize = 64 << 20;

#[derive(Default)]
struct Chunk {
    names: Vec<String>,
    held: Vec<(usize, Vec<u8>)>,
}

type Chunks = std::sync::mpsc::SyncSender<Result<Chunk>>;

pub fn call_over<F>(assembly: &str, min_length: usize, mut batch: F) -> Result<Vec<String>>
where
    F: FnMut(Vec<Orf>) -> Result<()>,
{
    let predictor =
        MetaPredictor::with_config_and_thread_pool(ProdigalConfig::default(), pool::get())
            .map_err(|error| anyhow!("gene finder: {error:?}"))?;

    let (sender, receiver) = sync_channel::<Result<Chunk>>(1);
    let held_assembly = assembly.to_string();
    let reader = thread::spawn(move || {
        if let Err(error) = read_chunks(&held_assembly, min_length, min_length, &sender) {
            let _ = sender.send(Err(error));
        }
    });

    let progress = crate::progress::spinning(crate::progress::Stage::CallingGenes);
    let mut names = Vec::new();
    let mut outcome = Ok(());
    for chunk in receiver {
        let mut chunk = match chunk {
            Ok(chunk) => chunk,
            Err(error) => {
                outcome = Err(error);
                break;
            }
        };
        names.append(&mut chunk.names);
        progress.set_message(format!("{} contigs", names.len()));
        if chunk.held.is_empty() {
            continue;
        }
        outcome = call(&predictor, &chunk.held).and_then(&mut batch);
        if outcome.is_err() {
            break;
        }
    }
    progress.finish_and_clear();
    reader
        .join()
        .map_err(|_| anyhow!("the assembly reader panicked"))?;
    outcome.map(|()| names)
}

fn read_chunks(
    assembly: &str,
    min_length: usize,
    called_from: usize,
    sender: &Chunks,
) -> Result<()> {
    let mut reader = parse_fastx_file(assembly)?;
    let mut chunk = Chunk::default();
    let mut seen = 0usize;
    let mut bases = 0usize;
    while let Some(record) = reader.next() {
        let record = record?;
        let sequence = record.seq();
        if sequence.len() < min_length {
            continue;
        }
        let name = crate::contig_id(record.id())?;
        if sequence.len() >= called_from {
            bases += sequence.len();
            chunk.held.push((seen, sequence.to_vec()));
        }
        chunk.names.push(name.to_string());
        seen += 1;
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

// The batch runs to the slowest contig, so the long ones have to start first.
fn call(predictor: &MetaPredictor, contigs: &[(usize, Vec<u8>)]) -> Result<Vec<Orf>> {
    let mut order = (0..contigs.len()).collect::<Vec<_>>();
    order.sort_unstable_by_key(|index| std::cmp::Reverse(contigs[*index].1.len()));
    let sequences = order
        .iter()
        .map(|index| contigs[*index].1.as_slice())
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
        .map(|((contig, sequence), genes)| {
            genes
                .into_iter()
                .map(|gene| {
                    let coding = &sequence[gene.begin - 1..gene.end];
                    let protein = match gene.strand {
                        Strand::Reverse => translate_reverse(coding, !gene.partial.1),
                        _ => translate(coding, !gene.partial.0),
                    };
                    Orf {
                        contig: *contig,
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
