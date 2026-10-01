use std::{
    collections::{HashMap, HashSet},
    io::{BufRead, Read},
    path::Path,
};

use anyhow::{Result, anyhow, bail};
use log::debug;
use ndarray::Array2;
use needletail::Sequence;
use rayon::prelude::*;

pub const DEFAULT_KMER_SIZE: usize = 4;

pub fn count_kmers(
    assembly: &str,
    output_directory: &str,
    n_contigs: usize,
    min_length: usize,
    kmer_size: usize,
    keep: bool,
) -> Result<KmerFrequencyTable> {
    KmerCounter {
        assembly: assembly.to_string(),
        output_directory: output_directory.to_string(),
        kmer_size,
        n_contigs,
        min_length,
        keep,
    }
    .run()
}

struct Block {
    kmer_size: usize,
    columns: Vec<u32>,
    width: usize,
}

struct KmerCounter {
    assembly: String,
    output_directory: String,
    kmer_size: usize,
    n_contigs: usize,
    min_length: usize,
    keep: bool,
}

impl KmerCounter {
    fn run(&mut self) -> Result<KmerFrequencyTable> {
        let output_file = Path::new(&self.output_directory)
            .join(format!("kmer_frequencies.k{}.tsv", self.kmer_size));
        if output_file.exists() {
            return KmerFrequencyTable::read(&output_file);
        }

        let block = block_of(self.kmer_size);
        let width = block.width;
        let mut kmer_table = Vec::with_capacity(self.n_contigs * width);
        let mut contig_names = Vec::with_capacity(self.n_contigs);
        let progress = crate::progress::spinning(crate::progress::Stage::CountingKmers);
        crate::kmers::measured(
            &self.assembly,
            self.min_length,
            |sequence| frequencies_of(sequence, &block),
            |chunk| {
                for (name, row) in chunk {
                    if let Some(row) = row {
                        contig_names.push(name);
                        kmer_table.extend(row);
                    }
                }
                progress.set_message(format!("{} contigs", contig_names.len()));
            },
        )?;
        progress.finish_and_clear();

        let kmer_array = Array2::from_shape_vec((contig_names.len(), width), kmer_table)?;

        let kmer_frequency_table =
            KmerFrequencyTable::new(self.kmer_size, kmer_array, contig_names);
        if self.keep {
            kmer_frequency_table.write(&output_file)?;
        }

        Ok(kmer_frequency_table)
    }
}

fn block_of(kmer_size: usize) -> Block {
    let canonical = canonical_index(kmer_size);
    Block {
        kmer_size,
        width: canonical.len(),
        columns: column_table(kmer_size, &canonical),
    }
}

pub fn halves(assembly: &str, names: &[&str], kmer_size: usize) -> Result<[Array2<f64>; 2]> {
    let block = block_of(kmer_size);
    let sequences = named_sequences(assembly, names)?;
    let split = names
        .iter()
        .map(|name| {
            let sequence = &sequences[*name];
            sequence.split_at(sequence.len() / 2)
        })
        .collect::<Vec<_>>();
    let [first, second] = [
        split.iter().map(|(first, _)| *first).collect::<Vec<_>>(),
        split.iter().map(|(_, second)| *second).collect::<Vec<_>>(),
    ];
    let lengths = |pieces: &[&[u8]]| pieces.iter().map(|piece| piece.len()).collect::<Vec<_>>();
    Ok([
        composition(&first, &lengths(&first), &block)?,
        composition(&second, &lengths(&second), &block)?,
    ])
}

/// The composition of the first `length` bases of each named contig, one row per piece.
pub fn prefixes(assembly: &str, pieces: &[(&str, usize)], kmer_size: usize) -> Result<Array2<f64>> {
    let block = block_of(kmer_size);
    let names = pieces.iter().map(|(name, _)| *name).collect::<Vec<_>>();
    let sequences = named_sequences(assembly, &names)?;
    let cut = pieces
        .iter()
        .map(|(name, length)| {
            let sequence = &sequences[*name];
            &sequence[..(*length).min(sequence.len())]
        })
        .collect::<Vec<_>>();
    let lengths = pieces.iter().map(|(_, length)| *length).collect::<Vec<_>>();
    composition(&cut, &lengths, &block)
}

// The wanted contigs are a band or a sample, never the assembly, so holding them to count in
// parallel costs little beside the rows they become.
fn named_sequences(assembly: &str, names: &[&str]) -> Result<HashMap<String, Vec<u8>>> {
    let wanted = names.iter().copied().collect::<HashSet<_>>();
    let mut found = HashMap::with_capacity(wanted.len());
    crate::kmers::each_named(
        assembly,
        |name| wanted.contains(name),
        |name, record| {
            found.insert(name.to_string(), record.normalize(false).into_owned());
            Ok(())
        },
    )?;
    if let Some(missing) = names.iter().find(|name| !found.contains_key(**name)) {
        bail!("{missing} is not in {assembly}");
    }
    Ok(found)
}

fn composition(pieces: &[&[u8]], lengths: &[usize], block: &Block) -> Result<Array2<f64>> {
    let rows = pieces
        .par_iter()
        .map(|piece| frequencies_of(piece, block))
        .collect::<Vec<_>>();
    let mut table = Array2::from_shape_vec((rows.len(), block.width), rows.concat())?;
    crate::kmers::clr::clr(&mut table, lengths, block.kmer_size)?;
    Ok(table)
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

// The forward code alone is enough because the column table already folds both strands, and
// a window with a base outside ACGT is skipped by counting only after `kmer_size` clean bases.
fn frequencies_of(sequence: &[u8], block: &Block) -> Vec<f64> {
    let mask = (1usize << (2 * block.kmer_size)) - 1;
    let mut counts = vec![0u32; block.width];
    let mut n_kmers = 0u32;
    let (mut code, mut clean) = (0usize, 0usize);
    for base in sequence {
        let Some(bits) = bits(*base) else {
            clean = 0;
            continue;
        };
        code = (code << 2 | bits) & mask;
        clean += 1;
        if clean >= block.kmer_size {
            counts[block.columns[code] as usize] += 1;
            n_kmers += 1;
        }
    }
    // An all-N contig holds no k-mer, and its zeros have to stay zeros rather than 0/0.
    let total = f64::from(n_kmers.max(1));
    counts
        .iter()
        .map(|count| f64::from(*count) / total)
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

fn bits(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn encode(kmer: &[u8]) -> Option<usize> {
    kmer.iter()
        .try_fold(0usize, |code, base| Some(code << 2 | bits(*base)?))
}

pub const KMER_SIZES: std::ops::RangeInclusive<i64> = 2..=6;

fn kmer_size_of(n_kmers: usize) -> Result<usize> {
    (*KMER_SIZES.start() as usize..=*KMER_SIZES.end() as usize)
        .find(|kmer_size| canonical_count(*kmer_size) == n_kmers)
        .ok_or_else(|| anyhow!("No k-mer size gives a table of {n_kmers} columns."))
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

    pub fn write<P: AsRef<Path>>(&self, output_file: P) -> Result<()> {
        crate::report_sink::write_atomically(output_file.as_ref(), |file| {
            let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_writer(file);
            for (contig_name, row) in self.contig_names.iter().zip(self.kmer_table.rows()) {
                writer.serialize((contig_name, row.into_iter().collect::<Vec<_>>()))?;
            }
            writer.flush()?;
            Ok(())
        })
    }

    /// Tables written before the switch to tabs are comma separated, and both still read.
    pub fn read<P: AsRef<Path>>(input_file: P) -> Result<Self> {
        let mut source = crate::get_file_reader(&input_file)?;
        let mut first = String::new();
        source.read_line(&mut first)?;
        let delimiter = if first.contains('\t') { b'\t' } else { b',' };
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(delimiter)
            .from_reader(first.as_bytes().chain(source));
        let mut contig_names = Vec::new();
        let mut kmer_table = Vec::new();
        for result in reader.deserialize() {
            let (contig_name, row): (String, Vec<f64>) = result?;
            contig_names.push(contig_name);
            kmer_table.push(row);
        }

        let Some(n_kmers) = kmer_table.first().map(Vec::len) else {
            bail!("{} holds no k-mer rows", input_file.as_ref().display());
        };
        debug!("Read n contigs {}", kmer_table.len());
        let kmer_size = kmer_size_of(n_kmers)?;

        let kmer_array = Array2::from_shape_vec(
            (contig_names.len(), n_kmers),
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
