use std::{
    collections::{HashMap, HashSet},
    path::Path,
};

use anyhow::{Result, bail};
use log::debug;
use ndarray::{Array, Array2, Axis};
use needletail::Sequence;
use rayon::prelude::*;

const DEFAULT_N_CONTIGS: usize = 10000;
pub const DEFAULT_KMER_SIZE: usize = 4;

/// Standard multiplicative replacement puts the substitute below the smallest observable
/// value rather than at it.
const REPLACEMENT_FRACTION: f64 = 0.65;

/// A contig cannot express a frequency below one count over its own kmer positions, so that
/// reciprocal is the floor a replacement has to sit under.
fn detection_limit(contig_length: usize, kmer_size: usize) -> f64 {
    let positions = contig_length.saturating_sub(kmer_size - 1).max(1);
    REPLACEMENT_FRACTION / positions as f64
}

fn median_replacement(contig_lengths: &[usize], kmer_size: usize) -> f64 {
    if contig_lengths.is_empty() {
        return f64::NAN;
    }
    let mut deltas = contig_lengths
        .iter()
        .map(|length| detection_limit(*length, kmer_size))
        .collect::<Vec<_>>();
    deltas.sort_by(f64::total_cmp);
    deltas[deltas.len() / 2]
}

pub fn count_kmers(
    assembly: &str,
    output_directory: &str,
    n_contigs: Option<usize>,
    kmer_size: usize,
) -> Result<KmerFrequencyTable> {
    KmerCounter::new(assembly, output_directory, n_contigs, kmer_size).run()
}

struct KmerCounter {
    assembly: String,
    output_directory: String,
    kmer_size: usize,
    n_contigs: Option<usize>,
}

impl KmerCounter {
    fn new(
        assembly: &str,
        output_directory: &str,
        n_contigs: Option<usize>,
        kmer_size: usize,
    ) -> Self {
        Self {
            assembly: assembly.to_string(),
            output_directory: output_directory.to_string(),
            kmer_size,
            n_contigs,
        }
    }

    fn run(&mut self) -> Result<KmerFrequencyTable> {
        let output_file = Path::new(&self.output_directory)
            .join(format!("kmer_frequencies.k{}.tsv", self.kmer_size));
        if output_file.exists() {
            return KmerFrequencyTable::read(&output_file);
        }

        let canonical_kmers = self.calculate_canonical_kmers();

        // use needletail to read in assembly and count canonical kmers
        // use ndarray to store kmer frequencies. 2D array with rows = contigs and columns = kmers
        let mut reader = needletail::parse_fastx_file(&self.assembly)?;

        let n_contigs = match self.n_contigs {
            Some(n) => n,
            None => DEFAULT_N_CONTIGS,
        };

        let mut kmer_table = Vec::with_capacity(n_contigs);
        let mut contig_names = Vec::with_capacity(n_contigs);
        let mut n_contigs = 0;
        while let Some(record) = reader.next() {
            let seqrec = record?;
            n_contigs += 1;
            let contig_name = std::str::from_utf8(seqrec.id())?.to_string();
            contig_names.push(contig_name);
            // normalize to make sure all the bases are consistently capitalized and
            // that we remove the newlines since this is FASTA
            let norm_seq = seqrec.normalize(false);
            // we make a reverse complemented copy of the sequence first for
            // `canonical_kmers` to draw the complemented sequences from.
            let rc = norm_seq.reverse_complement();
            // now we keep track of the number of AAAAs (or TTTTs via
            // canonicalization) in the file; note we also get the position (i.0;
            // in the event there were `N`-containing kmers that were skipped)
            // and whether the sequence was complemented (i.2) in addition to
            // the canonical kmer (i.1)
            let mut contig_kmer_counts = vec![0; canonical_kmers.len()];
            let mut n_kmers = 0;
            for (_, kmer, _) in norm_seq.canonical_kmers(self.kmer_size as u8, &rc) {
                // we need to calculate what the index of the kmer is in the
                // `contig_kmer_counts` vector; we do this by converting the
                // kmer to a base-4 number (A=0, C=1, G=2, T=3) and then
                // multiplying by 4^kmer_size-1, 4^kmer_size-2, etc. to get the
                // index
                let kmer_idx = if let Some(index) = canonical_kmers.get(kmer) {
                    *index
                } else {
                    // try the reverse complement?
                    let rc = kmer.reverse_complement();
                    if let Some(index) = canonical_kmers.get(&rc) {
                        *index
                    } else {
                        // we skip N-containing kmers
                        continue;
                    }
                };
                contig_kmer_counts[kmer_idx] += 1;
                n_kmers += 1;
            }
            // we need to convert the counts to frequencies
            let contig_kmer_freqs = contig_kmer_counts
                .iter()
                .map(|c| *c as f64 / n_kmers as f64)
                .collect::<Vec<f64>>();
            kmer_table.push(contig_kmer_freqs);
        }

        // convert kmer_table to Array2
        let kmer_array = Array2::from_shape_vec(
            (n_contigs, canonical_kmers.len()),
            kmer_table.into_iter().flatten().collect(),
        )?;

        let mut kmer_frequency_table = KmerFrequencyTable::new(
            self.kmer_size,
            kmer_array,
            contig_names,
            output_file.to_str().unwrap().to_string(),
        );
        kmer_frequency_table.write(&output_file)?;

        Ok(kmer_frequency_table)
    }

    fn calculate_canonical_kmers(&self) -> HashMap<Vec<u8>, usize> {
        canonical_index(self.kmer_size)
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

/// Canonical folding is not a power of four, so the width has to be matched against the class
/// count rather than inverted.
fn kmer_size_of(n_kmers: usize) -> Result<usize> {
    (1..=8)
        .find(|kmer_size| canonical_count(*kmer_size) == n_kmers)
        .ok_or_else(|| anyhow::anyhow!("No k-mer size gives a table of {n_kmers} columns."))
}

fn canonical_count(kmer_size: usize) -> usize {
    let palindromes = if kmer_size % 2 == 0 {
        4usize.pow(kmer_size as u32 / 2)
    } else {
        0
    };
    (4usize.pow(kmer_size as u32) + palindromes) / 2
}

/// increment a kmer to the next kmer in lexicographic order
fn increment_kmer(kmer: &mut [u8]) {
    // we start at the end of the kmer and increment the last base
    // if that base is a T, move the pointer to the next base and increment
    // that one, etc.
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
                    // we've reached the end of the kmer
                    break;
                } else {
                    // move to the next base
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
    pub(crate) table_path: String,
}

impl KmerFrequencyTable {
    pub fn new(
        kmer_size: usize,
        kmer_table: Array2<f64>,
        contig_names: Vec<String>,
        table_path: String,
    ) -> Self {
        Self {
            kmer_size,
            kmer_table,
            contig_names,
            table_path,
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
        // remove the contigs from the table
        let new_table = self
            .kmer_table
            .axis_iter(Axis(0))
            .enumerate()
            .filter_map(|(index, row)| {
                if indices_to_remove.contains(&index) {
                    None
                } else {
                    Some(row)
                }
            })
            .flat_map(|row| row.to_vec());
        let new_n_rows = self.kmer_table.nrows() - indices_to_remove.len();
        self.kmer_table = Array::from_iter(new_table)
            .into_shape_with_order((new_n_rows, self.kmer_table.ncols()))?;

        let filtered_contig_names = self
            .contig_names
            .iter()
            .enumerate()
            .filter_map(|(index, name)| {
                if indices_to_remove.contains(&index) {
                    Some(name.clone())
                } else {
                    None
                }
            })
            .collect::<HashSet<_>>();
        // remove the contigs from the contig names
        self.contig_names = self
            .contig_names
            .iter()
            .enumerate()
            .filter_map(|(index, name)| {
                if indices_to_remove.contains(&index) {
                    None
                } else {
                    Some(name.clone())
                }
            })
            .collect::<Vec<_>>();

        Ok(filtered_contig_names)
    }

    /// Write the kmer table to a file.
    pub fn write<P: AsRef<Path>>(&mut self, ouput_file: P) -> Result<()> {
        self.table_path = ouput_file.as_ref().to_str().unwrap().to_string();
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
            table_path: input_file.as_ref().to_str().unwrap().to_string(),
        })
    }

    /// Centre log ratio transform. Frequencies are compositional and the log cannot take the
    /// zeros short contigs leave, so they are replaced first. A single constant for the whole
    /// table is only right at about 18.5 kb, and 97% of a per-sample assembly sits below that.
    pub fn clr(&mut self, contig_lengths: &[usize]) -> Result<()> {
        let n_rows = self.kmer_table.nrows();
        let n_cols = self.kmer_table.ncols();
        let kmer_size = self.kmer_size;
        if contig_lengths.len() != n_rows {
            bail!(
                "Centre log ratio needs one length per row, got {} lengths for {} contigs.",
                contig_lengths.len(),
                n_rows
            );
        }

        let zeros: usize = self
            .kmer_table
            .iter()
            .filter(|value| **value <= 0.0)
            .count();
        debug!(
            "Zeros {:.4} of {} tetranucleotide cells, median replacement {:.3e}",
            zeros as f64 / (n_rows * n_cols) as f64,
            n_rows * n_cols,
            median_replacement(contig_lengths, kmer_size)
        );

        let new_array = (0..n_rows)
            .into_par_iter()
            .flat_map(|row_index| {
                let row = self.kmer_table.row(row_index);
                let row_sum = row.sum();
                let delta = detection_limit(contig_lengths[row_index], kmer_size);
                let n_zeros = row.iter().filter(|value| **value <= 0.0).count();
                let retained = (1.0 - n_zeros as f64 * delta).max(f64::MIN_POSITIVE);

                let replaced = (0..n_cols)
                    .map(|j| {
                        let value = row[[j]];
                        if value <= 0.0 {
                            delta
                        } else if row_sum > 0.0 {
                            value / row_sum * retained
                        } else {
                            delta
                        }
                    })
                    .collect::<Vec<_>>();

                let log_sum = replaced.iter().map(|value| value.ln()).sum::<f64>();
                let log_geometric_mean = log_sum / n_cols as f64;

                replaced
                    .into_iter()
                    .map(|value| value.ln() - log_geometric_mean)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        self.kmer_table = Array::from_shape_vec((n_rows, n_cols), new_array)?;
        Ok(())
    }
}
