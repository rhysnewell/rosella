use std::sync::Arc;

use anyhow::Result;
use needletail::parse_fastx_file;
use prodigal_rs::api::{META_PREDICTOR_STACK_SIZE, MetaPredictor, ProdigalConfig, Strand};
use rayon::prelude::*;

pub struct Orf {
    pub contig: usize,
    pub partial: bool,
    pub protein: String,
}

/// Contigs are called in chunks so the assembly is never held whole beside the proteins it
/// produces, which was the run's memory peak on a multi-sample assembly.
const CHUNK_BASES: usize = 64 << 20;

#[derive(Clone, Copy, Debug, Default)]
pub struct GeneRules {
    pub min_length: usize,
    pub model_depth: usize,
}

pub fn call_over<F>(
    assembly: &str,
    min_length: usize,
    rules: GeneRules,
    threads: usize,
    mut batch: F,
) -> Result<Vec<String>>
where
    F: FnMut(Vec<Orf>) -> Result<()>,
{
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(META_PREDICTOR_STACK_SIZE)
        .build()?;
    let config = ProdigalConfig {
        model_depth: rules.model_depth,
        ..ProdigalConfig::default()
    };
    let predictor = MetaPredictor::with_config_and_thread_pool(config, Arc::new(pool))
        .map_err(|error| anyhow!("gene finder: {error:?}"))?;
    let called_from = rules.min_length.max(min_length);

    let mut reader = parse_fastx_file(assembly)?;
    let mut names = Vec::new();
    let mut held: Vec<(usize, Vec<u8>)> = Vec::new();
    let mut bases = 0usize;
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
        if sequence.len() >= called_from {
            bases += sequence.len();
            held.push((names.len(), sequence.to_vec()));
        }
        names.push(name.to_string());
        if bases >= CHUNK_BASES {
            batch(call(&predictor, &held)?)?;
            held.clear();
            bases = 0;
        }
    }
    if !held.is_empty() {
        batch(call(&predictor, &held)?)?;
    }
    Ok(names)
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

pub fn translate(coding: &[u8], complete_start: bool) -> String {
    let mut protein = String::with_capacity(coding.len() / 3);
    for codon in coding.chunks_exact(3) {
        let residue = match (base(codon[0]), base(codon[1]), base(codon[2])) {
            (Some(a), Some(b), Some(c)) => CODONS[a * 16 + b * 4 + c] as char,
            _ => 'X',
        };
        protein.push(residue);
    }
    finish(protein, complete_start)
}

pub fn translate_reverse(coding: &[u8], complete_start: bool) -> String {
    let mut protein = String::with_capacity(coding.len() / 3);
    let mut end = coding.len();
    while end >= 3 {
        let residue = match (
            complement(coding[end - 1]),
            complement(coding[end - 2]),
            complement(coding[end - 3]),
        ) {
            (Some(a), Some(b), Some(c)) => CODONS[a * 16 + b * 4 + c] as char,
            _ => 'X',
        };
        protein.push(residue);
        end -= 3;
    }
    finish(protein, complete_start)
}
