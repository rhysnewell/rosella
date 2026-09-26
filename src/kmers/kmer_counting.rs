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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KmerSizes(Vec<usize>);

impl KmerSizes {
    pub fn as_slice(&self) -> &[usize] {
        &self.0
    }

    fn label(&self) -> String {
        self.0
            .iter()
            .map(|kmer_size| kmer_size.to_string())
            .collect::<Vec<_>>()
            .join("-")
    }
}

impl From<Vec<usize>> for KmerSizes {
    fn from(mut sizes: Vec<usize>) -> Self {
        sizes.sort_unstable();
        sizes.dedup();
        Self(sizes)
    }
}

impl std::str::FromStr for KmerSizes {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let sizes = value
            .split(',')
            .map(|part| {
                let kmer_size: i64 = part
                    .trim()
                    .parse()
                    .map_err(|_| format!("`{part}` is not a k-mer size"))?;
                if !KMER_SIZES.contains(&kmer_size) {
                    return Err(format!(
                        "`{kmer_size}` is outside {} to {}",
                        KMER_SIZES.start(),
                        KMER_SIZES.end()
                    ));
                }
                Ok(kmer_size as usize)
            })
            .collect::<Result<Vec<_>, Self::Err>>()?;
        if sizes.is_empty() {
            return Err("--kmer-size was given nothing".to_string());
        }
        Ok(Self::from(sizes))
    }
}

pub fn count_kmers(
    assembly: &str,
    output_directory: &str,
    n_contigs: Option<usize>,
    kmer_sizes: &KmerSizes,
    keep: bool,
) -> Result<KmerFrequencyTable> {
    KmerCounter::new(assembly, output_directory, n_contigs, kmer_sizes, keep).run()
}

struct Block {
    kmer_size: usize,
    columns: Vec<u32>,
    width: usize,
}

struct KmerCounter {
    assembly: String,
    output_directory: String,
    kmer_sizes: KmerSizes,
    n_contigs: Option<usize>,
    keep: bool,
}

impl KmerCounter {
    fn new(
        assembly: &str,
        output_directory: &str,
        n_contigs: Option<usize>,
        kmer_sizes: &KmerSizes,
        keep: bool,
    ) -> Self {
        Self {
            assembly: assembly.to_string(),
            output_directory: output_directory.to_string(),
            kmer_sizes: kmer_sizes.clone(),
            n_contigs,
            keep,
        }
    }

    fn run(&mut self) -> Result<KmerFrequencyTable> {
        let output_file = Path::new(&self.output_directory)
            .join(format!("kmer_frequencies.k{}.tsv", self.kmer_sizes.label()));
        if output_file.exists() {
            return KmerFrequencyTable::read(&output_file);
        }

        let (blocks, width) = blocks_of(&self.kmer_sizes);
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
                .map(|(_, sequence)| frequencies_of(sequence, &blocks, width))
                .collect::<Vec<_>>();
            for ((name, _), row) in chunk.iter().zip(frequencies) {
                contig_names.push(name.clone());
                kmer_table.extend(row);
            }
        }

        progress.finish_and_clear();

        let kmer_array = Array2::from_shape_vec((n_contigs, width), kmer_table)?;

        let mut kmer_frequency_table =
            KmerFrequencyTable::new(self.kmer_sizes.clone(), kmer_array, contig_names);
        if self.keep {
            kmer_frequency_table.write(&output_file)?;
        }

        Ok(kmer_frequency_table)
    }
}

fn blocks_of(kmer_sizes: &KmerSizes) -> (Vec<Block>, usize) {
    let blocks = kmer_sizes
        .as_slice()
        .iter()
        .map(|kmer_size| {
            let canonical = canonical_index(*kmer_size);
            Block {
                kmer_size: *kmer_size,
                width: canonical.len(),
                columns: column_table(*kmer_size, &canonical),
            }
        })
        .collect::<Vec<_>>();
    let width = blocks.iter().map(|block| block.width).sum::<usize>();
    (blocks, width)
}

pub fn halves(assembly: &str, names: &[&str], kmer_sizes: &KmerSizes) -> Result<[Array2<f64>; 2]> {
    let (blocks, width) = blocks_of(kmer_sizes);
    let wanted = names
        .iter()
        .enumerate()
        .map(|(at, name)| (*name, at))
        .collect::<HashMap<_, _>>();
    let mut rows = [vec![Vec::new(); names.len()], vec![Vec::new(); names.len()]];
    let mut lengths = [vec![0; names.len()], vec![0; names.len()]];
    let mut reader = needletail::parse_fastx_file(assembly)?;
    while let Some(record) = reader.next() {
        let record = record?;
        let Some(at) = wanted.get(crate::contig_id(record.id())?) else {
            continue;
        };
        let sequence = record.normalize(false);
        let (first, second) = sequence.split_at(sequence.len() / 2);
        for (side, half) in [first, second].into_iter().enumerate() {
            rows[side][*at] = frequencies_of(half, &blocks, width);
            lengths[side][*at] = half.len();
        }
    }
    if let Some(at) = rows[0].iter().position(Vec::is_empty) {
        bail!("{} is not in {assembly}", names[at]);
    }
    let [first, second] = rows.map(|held| {
        Array2::from_shape_vec((names.len(), width), held.concat())
            .expect("every row is one block set wide")
    });
    let mut tables = [first, second];
    for (table, lengths) in tables.iter_mut().zip(&lengths) {
        crate::kmers::clr::clr(table, lengths, kmer_sizes.as_slice())?;
    }
    Ok(tables)
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

fn frequencies_of(sequence: &[u8], blocks: &[Block], width: usize) -> Vec<f64> {
    let reverse = sequence.reverse_complement();
    let mut row = Vec::with_capacity(width);
    for block in blocks {
        let mut counts = vec![0u32; block.width];
        let mut n_kmers = 0u32;
        for (_, kmer, _) in sequence.canonical_kmers(block.kmer_size as u8, &reverse) {
            let Some(code) = encode(kmer) else { continue };
            counts[block.columns[code] as usize] += 1;
            n_kmers += 1;
        }
        row.extend(counts.iter().map(|count| *count as f64 / n_kmers as f64));
    }
    row
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

/// Each width is larger than every smaller width summed, so taking the largest that fits is
/// the only decomposition and the concatenation is recovered without carrying the list.
fn kmer_sizes_of(n_kmers: usize) -> Result<KmerSizes> {
    let mut left = n_kmers;
    let mut sizes = Vec::new();
    for kmer_size in (*KMER_SIZES.start() as usize..=*KMER_SIZES.end() as usize).rev() {
        let width = canonical_count(kmer_size);
        if width <= left {
            left -= width;
            sizes.push(kmer_size);
        }
    }
    if left != 0 {
        bail!("No set of k-mer sizes gives a table of {n_kmers} columns.");
    }
    Ok(KmerSizes::from(sizes))
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
    pub(crate) kmer_sizes: KmerSizes,
    pub kmer_table: Array2<f64>,
    pub(crate) contig_names: Vec<String>,
}

impl KmerFrequencyTable {
    pub fn new(kmer_sizes: KmerSizes, kmer_table: Array2<f64>, contig_names: Vec<String>) -> Self {
        Self {
            kmer_sizes,
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
            .from_reader(crate::get_file_reader(&input_file)?);
        let mut contig_names = Vec::new();
        let mut kmer_table = Vec::new();
        for result in reader.deserialize() {
            let (contig_name, row): (String, Vec<f64>) = result?;
            contig_names.push(contig_name);
            kmer_table.push(row);
        }

        let n_kmers = kmer_table[0].len();
        debug!("Read n contigs {}", kmer_table.len());
        let kmer_sizes = kmer_sizes_of(n_kmers)?;

        let kmer_array = Array2::from_shape_vec(
            (contig_names.len(), kmer_table[0].len()),
            kmer_table.into_iter().flatten().collect(),
        )?;

        Ok(Self {
            kmer_sizes,
            kmer_table: kmer_array,
            contig_names,
        })
    }

    pub fn kmer_sizes(&self) -> KmerSizes {
        self.kmer_sizes.clone()
    }

    pub fn clr(&mut self, contig_lengths: &[usize]) -> Result<()> {
        crate::kmers::clr::clr(
            &mut self.kmer_table,
            contig_lengths,
            self.kmer_sizes.as_slice(),
        )
    }
}
