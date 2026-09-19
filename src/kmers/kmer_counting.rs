use std::{
    collections::{HashMap, HashSet},
    path::Path,
};

use anyhow::{Result, bail};
use log::debug;
use ndarray::Array2;
use needletail::Sequence;
use rayon::prelude::*;

const DEFAULT_N_CONTIGS: usize = 10000;
pub const DEFAULT_KMER_SIZE: usize = 4;

pub fn count_kmers(
    assembly: &str,
    output_directory: &str,
    n_contigs: Option<usize>,
    kmer_size: usize,
    keep: bool,
) -> Result<KmerFrequencyTable> {
    KmerCounter::new(assembly, output_directory, n_contigs, kmer_size, keep).run()
}

struct KmerCounter {
    assembly: String,
    output_directory: String,
    kmer_size: usize,
    n_contigs: Option<usize>,
    keep: bool,
}

impl KmerCounter {
    fn new(
        assembly: &str,
        output_directory: &str,
        n_contigs: Option<usize>,
        kmer_size: usize,
        keep: bool,
    ) -> Self {
        Self {
            assembly: assembly.to_string(),
            output_directory: output_directory.to_string(),
            kmer_size,
            n_contigs,
            keep,
        }
    }

    fn run(&mut self) -> Result<KmerFrequencyTable> {
        let output_file = Path::new(&self.output_directory)
            .join(format!("kmer_frequencies.k{}.tsv", self.kmer_size));
        if output_file.exists() {
            return KmerFrequencyTable::read(&output_file);
        }

        let canonical = canonical_index(self.kmer_size);
        let width = canonical.len();
        let columns = column_table(self.kmer_size, &canonical);
        let mut reader = needletail::parse_fastx_file(&self.assembly)?;

        let expected = self.n_contigs.unwrap_or(DEFAULT_N_CONTIGS);
        let mut kmer_table = Vec::with_capacity(expected * width);
        let mut contig_names = Vec::with_capacity(expected);
        let mut chunk: Vec<(String, Vec<u8>)> = Vec::with_capacity(CHUNK);
        let mut n_contigs = 0;
        let progress = crate::progress::spinning(crate::progress::Stage::CountingKmers);
        loop {
            chunk.clear();
            while chunk.len() < CHUNK {
                let Some(record) = reader.next() else { break };
                let seqrec = record?;
                let name = crate::contig_id(seqrec.id())?.to_string();
                chunk.push((name, seqrec.normalize(false).into_owned()));
            }
            if chunk.is_empty() {
                break;
            }
            n_contigs += chunk.len();
            progress.set_message(format!("{n_contigs} contigs"));
            let frequencies = chunk
                .par_iter()
                .map(|(_, sequence)| frequencies_of(sequence, self.kmer_size, &columns, width))
                .collect::<Vec<_>>();
            for ((name, _), row) in chunk.iter().zip(frequencies) {
                contig_names.push(name.clone());
                kmer_table.extend(row);
            }
        }

        progress.finish_and_clear();

        let kmer_array = Array2::from_shape_vec((n_contigs, width), kmer_table)?;

        let mut kmer_frequency_table =
            KmerFrequencyTable::new(self.kmer_size, kmer_array, contig_names);
        if self.keep {
            kmer_frequency_table.write(&output_file)?;
        }

        Ok(kmer_frequency_table)
    }
}

/// Every k-mer folded onto the lexicographically smaller of itself and its reverse complement,
/// mapped to its sorted position, which is the table's column order.
pub fn canonical_index(kmer_size: usize) -> HashMap<Vec<u8>, usize> {
    let mut canonical_kmers = HashSet::with_capacity(4usize.pow(kmer_size as u32));

    let mut kmer = vec![b'A'; kmer_size];
    for _ in 0..4usize.pow(kmer_size as u32) {
        let rc_kmer = kmer.reverse_complement();
        if !canonical_kmers.contains(&kmer) && !canonical_kmers.contains(&rc_kmer) {
            if kmer < rc_kmer {
                canonical_kmers.insert(kmer.clone());
            } else {
                canonical_kmers.insert(rc_kmer.clone());
            }
        }

        increment_kmer(&mut kmer);
    }

    let mut canonical_kmers = canonical_kmers.into_iter().collect::<Vec<_>>();
    canonical_kmers.par_sort_unstable();
    canonical_kmers
        .into_par_iter()
        .enumerate()
        .map(|(i, k)| (k, i))
        .collect()
}

/// Contigs are read in chunks rather than whole so the parallel count does not hold the
/// assembly in memory beside the table it is filling.
const CHUNK: usize = 512;

fn frequencies_of(sequence: &[u8], kmer_size: usize, columns: &[u32], width: usize) -> Vec<f64> {
    let reverse = sequence.reverse_complement();
    let mut counts = vec![0u32; width];
    let mut n_kmers = 0u32;
    for (_, kmer, _) in sequence.canonical_kmers(kmer_size as u8, &reverse) {
        let Some(code) = encode(kmer) else { continue };
        counts[columns[code] as usize] += 1;
        n_kmers += 1;
    }
    counts
        .iter()
        .map(|count| *count as f64 / n_kmers as f64)
        .collect()
}

/// Every two-bit encoding indexed straight to its canonical column, so counting costs no hash
/// per base and needs no reverse complement lookup for the half that folds.
fn column_table(kmer_size: usize, canonical: &HashMap<Vec<u8>, usize>) -> Vec<u32> {
    let mut table = vec![0u32; 4usize.pow(kmer_size as u32)];
    let mut kmer = vec![b'A'; kmer_size];
    for _ in 0..table.len() {
        let column = canonical
            .get(&kmer)
            .or_else(|| canonical.get(&kmer.reverse_complement()))
            .expect("canonical folding covers every kmer");
        table[encode(&kmer).expect("generated kmers hold no ambiguity")] = *column as u32;
        increment_kmer(&mut kmer);
    }
    table
}

fn encode(kmer: &[u8]) -> Option<usize> {
    kmer.iter().try_fold(0usize, |code, base| {
        let bits = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        Some(code << 2 | bits)
    })
}

pub const KMER_SIZES: std::ops::RangeInclusive<i64> = 2..=6;

fn kmer_size_of(n_kmers: usize) -> Result<usize> {
    for kmer_size in *KMER_SIZES.start() as usize..=*KMER_SIZES.end() as usize {
        if canonical_count(kmer_size) == n_kmers {
            return Ok(kmer_size);
        }
    }
    bail!("No k-mer size gives a table of {n_kmers} columns.")
}

pub fn canonical_count(kmer_size: usize) -> usize {
    let palindromes = if kmer_size % 2 == 0 {
        4usize.pow(kmer_size as u32 / 2)
    } else {
        0
    };
    (4usize.pow(kmer_size as u32) + palindromes) / 2
}

fn increment_kmer(kmer: &mut [u8]) {
    let mut i = kmer.len() - 1;
    loop {
        match kmer[i] {
            b'A' => {
                kmer[i] = b'C';
                break;
            }
            b'C' => {
                kmer[i] = b'G';
                break;
            }
            b'G' => {
                kmer[i] = b'T';
                break;
            }
            b'T' => {
                kmer[i] = b'A';
                if i == 0 {
                    break;
                } else {
                    i -= 1;
                }
            }
            _ => unreachable!(),
        }
    }
}

pub struct KmerFrequencyTable {
    pub(crate) kmer_size: usize,
    pub kmer_table: Array2<f64>,
    pub(crate) contig_names: Vec<String>,
}

impl KmerFrequencyTable {
    pub fn new(kmer_size: usize, kmer_table: Array2<f64>, contig_names: Vec<String>) -> Self {
        Self {
            kmer_size,
            kmer_table,
            contig_names,
        }
    }

    pub fn filter_by_name(&mut self, to_filter: &HashSet<String>) -> Result<HashSet<String>> {
        // find the indices of the contigs that are too small
        let indices_to_remove = self
            .contig_names
            .iter()
            .enumerate()
            .filter_map(|(index, name)| {
                if to_filter.contains(name) {
                    Some(index)
                } else {
                    None
                }
            })
            .collect::<HashSet<_>>();

        self.filter_by_index(&indices_to_remove)
    }

    pub fn filter_by_index(
        &mut self,
        indices_to_remove: &HashSet<usize>,
    ) -> Result<HashSet<String>> {
        let removed = crate::rows::dropped_names(&self.contig_names, indices_to_remove);
        self.kmer_table = crate::rows::keep_rows(&self.kmer_table, indices_to_remove)?;
        self.contig_names = crate::rows::keep(&self.contig_names, indices_to_remove);
        Ok(removed)
    }

    /// Write the kmer table to a file.
    pub fn write<P: AsRef<Path>>(&mut self, ouput_file: P) -> Result<()> {
        let mut writer = csv::Writer::from_path(ouput_file)?;
        // we won't write a header for this file.
        for (contig_name, row) in self.contig_names.iter().zip(self.kmer_table.rows()) {
            writer.serialize((contig_name, row.into_iter().collect::<Vec<_>>()))?;
        }
        writer.flush()?;

        Ok(())
    }

    /// Read a kmer table from a file.
    pub fn read<P: AsRef<Path>>(input_file: P) -> Result<Self> {
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_path(&input_file)?;
        let mut contig_names = Vec::new();
        let mut kmer_table = Vec::new();
        for result in reader.deserialize() {
            let (contig_name, row): (String, Vec<f64>) = result?;
            contig_names.push(contig_name);
            kmer_table.push(row);
        }

        let n_kmers = kmer_table[0].len();
        debug!("Read n contigs {}", kmer_table.len());
        let kmer_size = kmer_size_of(n_kmers)?;

        let kmer_array = Array2::from_shape_vec(
            (contig_names.len(), kmer_table[0].len()),
            kmer_table.into_iter().flatten().collect(),
        )?;

        Ok(Self {
            kmer_size,
            kmer_table: kmer_array,
            contig_names,
        })
    }

    pub fn kmer_size(&self) -> usize {
        self.kmer_size
    }

    pub fn clr(&mut self, contig_lengths: &[usize]) -> Result<()> {
        crate::kmers::clr::clr(&mut self.kmer_table, contig_lengths, self.kmer_size)
    }
}
