use std::path::Path;

use anyhow::Result;
use ndarray::Array2;
use needletail::Sequence;

const DEFAULT_N_CONTIGS: usize = 10000;

pub fn count_kmers(m: &clap::ArgMatches, n_contigs: Option<usize>) -> Result<KmerFrequencyTable> {
    let mut kmer_counter = KmerCounter::new(m, n_contigs)?;
    kmer_counter.run()
}

struct KmerCounter {
    assembly: String,
    output_directory: String,
    kmer_size: usize,
    n_contigs: Option<usize>,
}

impl KmerCounter {
    pub fn new(m: &clap::ArgMatches, n_contigs: Option<usize>) -> Result<Self> {
        let assembly = m.get_one::<String>("assembly").unwrap().clone();
        let output_directory = m.get_one::<String>("output-directory").unwrap().clone();
        let kmer_size = *m.get_one::<usize>("kmer-size").unwrap();

        Ok(
            Self {
                assembly,
                output_directory,
                kmer_size,
                n_contigs,
            }
        )
    }

    pub fn run(&mut self) -> Result<KmerFrequencyTable> {
        let output_file = Path::new(&self.output_directory).join("kmer_frequencies.tsv");
        // check if output file exists
        if output_file.exists() {
            // if it does, read it in
            let kmer_table = KmerFrequencyTable::read(&output_file)?;
            return Ok(kmer_table)
        }

        // use needletail to read in assembly and count canonical kmers
        // use ndarray to store kmer frequencies. 2D array with rows = contigs and columns = kmers
        let mut reader = needletail::parse_fastx_file(&self.assembly)?;
        
        let n_contigs = match self.n_contigs {
            Some(n) => n,
            None => DEFAULT_N_CONTIGS,
        };

        let mut kmer_table = Vec::with_capacity(n_contigs);
        let mut contig_names = Vec::with_capacity(n_contigs);
        while let Some(record) = reader.next() {
            let seqrec = record?;
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
            let mut contig_kmer_counts = vec![0; 4usize.pow(self.kmer_size as u32)];
            let mut n_kmers = 0;
            for (_, kmer, _) in norm_seq.canonical_kmers(self.kmer_size as u8, &rc) {
                // we need to calculate what the index of the kmer is in the
                // `contig_kmer_counts` vector; we do this by converting the
                // kmer to a base-4 number (A=0, C=1, G=2, T=3) and then
                // multiplying by 4^kmer_size-1, 4^kmer_size-2, etc. to get the
                // index
                let kmer_idx = kmer
                    .iter()
                    .map(|b| match b {
                        b'A' | b'a' => 0,
                        b'C' | b'c' => 1,
                        b'G' | b'g' => 2,
                        b'T' | b't' => 3,
                        _ => unreachable!(),
                    })
                    .enumerate()
                    .map(|(i, b)| b * 4usize.pow((self.kmer_size - 1 - i) as u32))
                    .sum::<usize>();
                contig_kmer_counts[kmer_idx] += 1;
                n_kmers += 1;
            };
            // we need to convert the counts to frequencies
            let contig_kmer_freqs = contig_kmer_counts
                .iter()
                .map(|c| *c as f32 / n_kmers as f32)
                .collect::<Vec<f32>>();
            kmer_table.push(contig_kmer_freqs);
        };

        // convert kmer_table to Array2
        let kmer_array = Array2::from_shape_vec(
            (n_contigs, 4usize.pow(self.kmer_size as u32)), 
            kmer_table.into_iter().flatten().collect())?;
        
        let kmer_frequency_table = KmerFrequencyTable::new(self.kmer_size, kmer_array, contig_names);
        kmer_frequency_table.write(&output_file)?;

        Ok(kmer_frequency_table)
    }
}

pub struct KmerFrequencyTable {
    kmer_size: usize,
    kmer_table: Array2<f32>,
    contig_names: Vec<String>,
}

impl KmerFrequencyTable {
    pub fn new(kmer_size: usize, kmer_table: Array2<f32>, contig_names: Vec<String>) -> Self {
        Self {
            kmer_size,
            kmer_table,
            contig_names,
        }
    }

    /// Write the kmer table to a file.
    pub fn write<P: AsRef<Path>>(&self, ouput_file: P) -> Result<()> {
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
        let mut reader = csv::Reader::from_path(input_file)?;
        let mut contig_names = Vec::new();
        let mut kmer_table = Vec::new();
        for result in reader.deserialize() {
            let (contig_name, row): (String, Vec<f32>) = result?;
            contig_names.push(contig_name);
            kmer_table.push(row);
        }

        let n_kmers = kmer_table[0].len();
        let kmer_size = (n_kmers as f32).log(4.0).round() as usize;

        let kmer_array = Array2::from_shape_vec(
            (contig_names.len(), kmer_table[0].len()), 
            kmer_table.into_iter().flatten().collect())?;

        Ok(
            Self {
                kmer_size,
                kmer_table: kmer_array,
                contig_names,
            }
        )
    }
}