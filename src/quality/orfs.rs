use std::sync::Arc;

use anyhow::Result;
use needletail::parse_fastx_file;
use prodigal_rs::api::{META_PREDICTOR_STACK_SIZE, MetaPredictor, ProdigalConfig, Strand};

pub struct Orf {
    pub contig: usize,
    pub partial: bool,
    pub protein: String,
}

type Contigs = (Vec<String>, Vec<(usize, Vec<u8>)>);

pub fn read_over(assembly: &str, min_length: usize) -> Result<Contigs> {
    let mut reader = parse_fastx_file(assembly)?;
    let mut names = Vec::new();
    let mut contigs = Vec::new();
    while let Some(record) = reader.next() {
        let record = record?;
        let sequence = record.seq();
        if sequence.len() < min_length {
            continue;
        }
        let name = std::str::from_utf8(record.id())?
            .split_whitespace()
            .next()
            .unwrap_or_default();
        contigs.push((names.len(), sequence.to_vec()));
        names.push(name.to_string());
    }
    Ok((names, contigs))
}

pub fn call(contigs: &[(usize, Vec<u8>)], threads: usize) -> Result<Vec<Orf>> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(META_PREDICTOR_STACK_SIZE)
        .build()?;
    let predictor =
        MetaPredictor::with_config_and_thread_pool(ProdigalConfig::default(), Arc::new(pool))
            .map_err(|error| anyhow!("gene finder: {error:?}"))?;
    let sequences = contigs
        .iter()
        .map(|(_, sequence)| sequence.as_slice())
        .collect::<Vec<_>>();
    let batches = predictor
        .predict_batch(&sequences)
        .map_err(|error| anyhow!("gene finder: {error:?}"))?;

    let mut orfs = Vec::new();
    for ((contig, sequence), genes) in contigs.iter().zip(batches) {
        for gene in genes {
            let coding = &sequence[gene.begin - 1..gene.end];
            let protein = match gene.strand {
                Strand::Reverse => translate(&reverse_complement(coding), !gene.partial.1),
                _ => translate(coding, !gene.partial.0),
            };
            orfs.push(Orf {
                contig: *contig,
                partial: gene.partial.0 || gene.partial.1,
                protein,
            });
        }
    }
    Ok(orfs)
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

pub fn translate(coding: &[u8], complete_start: bool) -> String {
    let mut protein = String::with_capacity(coding.len() / 3);
    for codon in coding.chunks_exact(3) {
        let residue = match (base(codon[0]), base(codon[1]), base(codon[2])) {
            (Some(a), Some(b), Some(c)) => CODONS[a * 16 + b * 4 + c] as char,
            _ => 'X',
        };
        protein.push(residue);
    }
    if protein.ends_with('*') {
        protein.pop();
    }
    if complete_start && !protein.is_empty() {
        protein.replace_range(0..1, "M");
    }
    protein
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|byte| match byte.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' | b'U' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => other,
        })
        .collect()
}
