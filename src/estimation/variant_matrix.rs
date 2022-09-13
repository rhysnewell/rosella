use crate::estimation::contig::Elem;
use crate::external_command_checker;
use bio::alignment::sparse::hash_kmers;
use bio::io::fasta::IndexedReader;
use bird_tool_utils::command;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use itertools::izip;
use itertools::Itertools;
use model::variants::*;
use needletail::{parse_fastx_file, FastxReader, Sequence};
use rayon::prelude::*;
use scoped_threadpool::Pool;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::str;
use std::sync::mpsc::channel;
use std::sync::{Arc, Mutex};
use utils::{fetch_contig_from_reference, generate_faidx, read_sequence_to_vec, ReadType, Samples};

#[derive(Debug, Clone)]
/// Container for all variants within a genome and associated clusters
pub enum VariantMatrix<'b> {
    VariantContigMatrix {
        samples: Samples,
        average_genotypes: HashMap<i32, Vec<f64>>,
        // TID, Position, Variant Type, Base
        all_variants: HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>,
        // Placeholder hashmap for the depths of each contig for a sample
        // Deleted after use
        target_names: BTreeMap<i32, String>,
        target_lengths: HashMap<i32, u64>,
        target_ids: HashMap<String, i32>,
        kfrequencies: BTreeMap<Vec<u8>, Vec<usize>>,
        kmerfrequencies: HashMap<i32, HashMap<Vec<u8>, usize>>,
        present_kmers: BTreeSet<Vec<u8>>,
        variant_counts: HashMap<i32, usize>,
        variant_sums: HashMap<i32, Vec<Vec<f64>>>,
        variant_info: Vec<Var>,
        geom_mean_var: Vec<f64>,
        geom_mean_dep: Vec<f64>,
        geom_mean_frq: Vec<f64>,
        variant_rates: HashMap<i32, HashMap<Var, Vec<f32>>>,
        final_bins: HashMap<usize, Vec<i32>>,
        _m: std::marker::PhantomData<&'b ()>,
        //        pred_variants_all: HashMap<usize, HashMap<i32, HashMap<i32, HashSet<String>>>>,
    },
}

impl<'b> VariantMatrix<'b> {
    pub fn new_matrix(
        short_sample_count: usize,
        long_sample_count: usize,
        references: &[String],
    ) -> VariantMatrix<'b> {

        let mut target_names = BTreeMap::new();
        let mut target_lengths = HashMap::new();
        let mut target_ids = HashMap::new();
        references.into_iter().for_each(|reference| {
            let reference = bio::io::fasta::Index::from_file(
                &format!("{}.fai", &reference)
            ).unwrap();

            for (tid, seq) in reference.sequences().into_iter().enumerate() {
                target_names.insert(tid as i32, seq.name.clone());
                target_lengths.insert(tid as i32, seq.len);
                target_ids.insert(seq.name, tid as i32);
            }
        });


        VariantMatrix::VariantContigMatrix {
            samples: Samples::new_with_counts(short_sample_count, long_sample_count),
            average_genotypes: HashMap::new(),
            all_variants: HashMap::new(),
            target_names,
            target_lengths,
            target_ids,
            kfrequencies: BTreeMap::new(),
            kmerfrequencies: HashMap::new(),
            present_kmers: BTreeSet::new(),
            variant_counts: HashMap::new(),
            variant_sums: HashMap::new(),
            variant_info: Vec::new(),
            geom_mean_var: Vec::new(),
            geom_mean_dep: Vec::new(),
            geom_mean_frq: Vec::new(),
            variant_rates: HashMap::new(),
            final_bins: HashMap::new(),
            _m: std::marker::PhantomData::default(),
        }
    }
}

pub trait VariantMatrixFunctions {
    fn setup(&mut self);

    /// returns sample name by index
    fn get_sample_name(&self, sample_idx: usize, read_type: ReadType) -> &str;

    /// returns the number of contigs
    fn get_n_contigs(&self) -> u32;

    /// returns the contig lengths
    fn get_contig_lengths(&self) -> HashMap<i32, u64>;

    fn add_sample_name(&mut self, sample_name: String, sample_idx: usize, read_type: ReadType);

    fn add_info(&mut self, tid: usize, target_name: Vec<u8>, target_len: u64);

    /// Takes [VariantStats](contig_variants/VariantStats) struct for single contig and adds to
    /// [VariantMatrix](VariantMatrix)
    fn add_contig(
        &mut self,
        tid: i32,
        coverage_info: Vec<f64>,
        sample_count: usize,
        sample_idx: usize,
        read_type: ReadType,
    );

    /// Calculates the kmer frequencies for a given contig
    fn calc_kmer_frequencies(
        &mut self,
        kmer_size: u8,
        references: &[String],
        contig_count: usize,
        min_contig_size: u64,
        pb: &Elem,
    );

    /// Calculates the kmer frequencies for a given contig
    fn calc_kmer_frequencies_single(
        &mut self,
        tid: u32,
        kmer_size: usize,
        reference_file: &mut IndexedReader<File>,
        contig_count: usize,
    );

    fn write_kmer_table(&self, output: &str, min_contig_size: u64);

    fn write_coverage(&self, output: &str, table_type: ReadType);

    fn merge_matrices(
        &mut self,
        tid: u32,
        other_matrix: VariantMatrix,
        assembly_sample_count: usize,
    );

    fn bin_contigs(&self, output: &str, m: &clap::ArgMatches);

    fn refine_bins(&self, output: &str, m: &clap::ArgMatches);

    fn finalize_bins(&mut self, output: &str, m: &clap::ArgMatches);

    fn write_bins(&self, output: &str, reference: &str, mode: &str);

    fn read_inputs(
        &mut self,
        coverages_path: Option<&str>,
        long_coverage_path: Option<&str>,
        kmer_path: Option<&str>,
    );
}

#[allow(unused)]
impl VariantMatrixFunctions for VariantMatrix<'_> {
    fn setup(&mut self) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut samples,
                ref mut average_genotypes,
                ref mut all_variants,
                ref mut target_names,
                ref mut target_lengths,
                ref mut kfrequencies,
                ..
            } => {
                *samples = Samples::new();
                *average_genotypes = HashMap::new();
                *all_variants = HashMap::new();
                *target_names = BTreeMap::new();
                *target_lengths = HashMap::new();
                *kfrequencies = BTreeMap::new();
            }
        }
    }

    fn get_sample_name(&self, sample_idx: usize, read_type: ReadType) -> &str {
        match self {
            VariantMatrix::VariantContigMatrix { samples, .. } => match read_type {
                ReadType::Short => {
                    if samples.short[sample_idx].to_string().contains(".tmp") {
                        return &samples.short[sample_idx][15..];
                    } else {
                        return &samples.short[sample_idx];
                    }
                }
                ReadType::Long => {
                    if samples.long[sample_idx].to_string().contains(".tmp") {
                        return &samples.long[sample_idx][15..];
                    } else {
                        return &samples.long[sample_idx];
                    }
                }
                _ => panic!("ReadType {:?} not supported", read_type),
            },
        }
    }

    fn get_n_contigs(&self) -> u32 {
        match self {
            VariantMatrix::VariantContigMatrix { target_names, .. } => {
                target_names.len().clone() as u32
            }
        }
    }

    fn get_contig_lengths(&self) -> HashMap<i32, u64> {
        match self {
            VariantMatrix::VariantContigMatrix { target_lengths, .. } => target_lengths.clone(),
        }
    }

    fn add_sample_name(&mut self, sample_name: String, sample_idx: usize, read_type: ReadType) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut samples, ..
            } => {
                debug!(
                    "adding sample {} at index {} for {:?} read type",
                    &sample_name, &sample_idx, &read_type
                );
                match read_type {
                    ReadType::Short => samples.short[sample_idx] = sample_name,
                    ReadType::Long => samples.long[sample_idx] = sample_name,
                    _ => panic!("ReadType not supported"),
                }
            }
        }
    }

    fn add_info(&mut self, tid: usize, target_name: Vec<u8>, target_len: u64) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut target_names,
                ref mut target_lengths,
                ..
            } => {
                if !target_names.contains_key(&(tid as i32)) {
                    target_names
                        .entry(tid as i32)
                        .or_insert(String::from_utf8(target_name).unwrap()); // add contig name

                    target_lengths.entry(tid as i32).or_insert(target_len); // add contig length
                }
            }
        }
    }

    fn add_contig(
        &mut self,
        tid: i32,
        coverage_info: Vec<f64>,
        sample_count: usize,
        sample_idx: usize,
        read_type: ReadType,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut samples,
                ref mut all_variants,
                target_names,
                ..
            } => match read_type {
                ReadType::Short => {
                    let var = samples
                        .short_variances
                        .entry(tid)
                        .or_insert(vec![0.0 as f64; sample_count]);
                    var[sample_idx] = coverage_info[2];
                    let cov = samples
                        .short_coverages
                        .entry(tid)
                        .or_insert(vec![0.0 as f64; sample_count]);
                    cov[sample_idx] = coverage_info[1];
                }
                ReadType::Long => {
                    let var = samples
                        .long_variances
                        .entry(tid)
                        .or_insert(vec![0.0 as f64; sample_count]);
                    var[sample_idx] = coverage_info[2];
                    let cov = samples
                        .long_coverages
                        .entry(tid)
                        .or_insert(vec![0.0 as f64; sample_count]);
                    cov[sample_idx] = coverage_info[1];
                }
                _ => panic!("ReadType {:?} not supported", read_type),
            },
        }
    }

    fn calc_kmer_frequencies(
        &mut self,
        kmer_size: u8,
        references: &[String],
        contig_count: usize,
        min_contig_size: u64,
        pb: &Elem,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut kfrequencies,
                ref mut kmerfrequencies,
                ref mut present_kmers,
                ..
            } => {
                let mut tid = 0;

                for reference in references {
                    let mut reader =
                        parse_fastx_file(reference).expect("invalid path/file for assembly");


                    while let Some(record) = reader.next() {
                        let seqrec = record.expect("Invalid record");
                        if seqrec.num_bases() >= min_contig_size as usize {
                            // normalize to make sure all the bases are consistently capitalized and
                            // that we remove the newlines since this is FASTA
                            let norm_seq = seqrec.normalize(true);
                            // we make a reverse complemented copy of the sequence first for
                            // `canonical_kmers` to draw the complemented sequences from.
                            let rc = norm_seq.reverse_complement();
                            // now we keep track of the number of AAAAs (or TTTTs via
                            // canonicalization) in the file; note we also get the position (i.0;
                            // in the event there were `N`-containing kmers that were skipped)
                            // and whether the sequence was complemented (i.2) in addition to
                            // the canonical kmer (i.1)
                            let mut kmer_count = HashMap::new();

                            for (_, kmer, _) in norm_seq.canonical_kmers(kmer_size, &rc) {
                                let mut acc = kmer_count.entry(kmer.to_vec()).or_insert(0);
                                *acc += 1;
                            }

                            if present_kmers.len() < 136 {
                                let current_kmers =
                                    kmer_count.keys().cloned().collect::<BTreeSet<Vec<u8>>>();
                                present_kmers.par_extend(current_kmers);
                            }

                            kmerfrequencies.insert(tid, kmer_count);
                        }

                        tid += 1;
                        if tid % 100 == 0 {
                            pb.progress_bar.inc(100);
                            pb.progress_bar
                                .set_message(format!("Contigs kmers analyzed..."));
                            let pos = pb.progress_bar.position();
                            let len = pb.progress_bar.length();
                            if pos >= len {
                                pb.progress_bar
                                    .finish_with_message(format!("All contigs analyzed {}", "✔",));
                            }
                        }
                    }
                }
                pb.progress_bar
                    .finish_with_message(format!("All contigs analyzed {}", "✔",));
            }
        }
    }

    fn calc_kmer_frequencies_single(
        &mut self,
        tid: u32,
        kmer_size: usize,
        reference_file: &mut IndexedReader<File>,
        contig_count: usize,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut kmerfrequencies,
                ref mut present_kmers,
                target_names,
                target_lengths,
                ..
            } => {
                let mut ref_seq = Vec::new();
                let placeholder_name = "NOT_FOUND".to_string();

                let target_name = match target_names.get(&(tid as i32)) {
                    Some(name) => name,
                    None => &placeholder_name,
                };
                // Update all contig information
                fetch_contig_from_reference(reference_file, &target_name.as_bytes().to_vec());
                read_sequence_to_vec(
                    &mut ref_seq,
                    reference_file,
                    &target_name.as_bytes().to_vec(),
                );
                ref_seq.make_ascii_uppercase();
                let kmers = hash_kmers(&ref_seq[..], kmer_size);
                let dna_alphabet = bio::alphabets::dna::alphabet();
                // Get kmer counts in a contig
                let kmer_count = kmers
                    .into_par_iter()
                    .filter(|(kmer, _)| dna_alphabet.is_word(*kmer))
                    .map(|(kmer, _)| {
                        let mut acc = HashMap::new();
                        // Get canonical kmer, which ever of the current kmer and its RC
                        // come first alphabetically/lexographically
                        let rc = bio::alphabets::dna::revcomp(kmer);
                        let mut k_list = vec![kmer, &rc[..]];
                        k_list.sort();

                        let kmer = k_list[0];

                        let mut k = acc.entry(kmer.to_vec()).or_insert(0);
                        *k += 1;
                        acc
                    })
                    .reduce(
                        || HashMap::new(),
                        |m1, m2| {
                            m2.iter().fold(m1, |mut acc, (k, vs)| {
                                acc.entry(k.clone()).or_insert(*vs);
                                acc
                            })
                        },
                    );

                if present_kmers.len() < 136 {
                    let current_kmers = kmer_count.keys().cloned().collect::<BTreeSet<Vec<u8>>>();
                    present_kmers.par_extend(current_kmers);
                }

                kmerfrequencies.insert(tid as i32, kmer_count);
            }
        }
    }

    fn write_kmer_table(&self, output: &str, min_contig_size: u64) {
        match self {
            VariantMatrix::VariantContigMatrix {
                kmerfrequencies,
                present_kmers,
                target_names,
                target_lengths,
                ..
            } => {
                let file_name = format!("{}/rosella_kmer_table.tsv", &output,);

                let file_path = Path::new(&file_name);

                let mut file_open = match File::create(file_path) {
                    Ok(file) => file,
                    Err(e) => {
                        println!("Cannot create file {:?}", e);
                        std::process::exit(1)
                    }
                };

                // Write kmers in headers
                write!(file_open, "{}", "contigName").unwrap();
                write!(file_open, "\t{}", "contigLen").unwrap();
                for kmer in present_kmers.iter() {
                    write!(file_open, "\t{}", str::from_utf8(&kmer[..]).unwrap()).unwrap();
                }

                write!(file_open, "\n").unwrap();

                for (idx, contig) in target_names.iter() {
                    // Write the length
                    let length = target_lengths.get(idx).unwrap();

                    if length >= &min_contig_size {
                        write!(file_open, "{}", contig).unwrap();

                        write!(file_open, "\t{}", length).unwrap();

                        let counts = kmerfrequencies.get(idx).unwrap();
                        for kmer in present_kmers.iter() {
                            let count = counts.get(kmer).unwrap_or(&0);
                            write!(file_open, "\t{}", count).unwrap();
                        }
                        write!(file_open, "\n").unwrap();
                    }
                }
            }
        }
    }

    fn write_coverage(&self, output: &str, table_type: ReadType) {
        match self {
            VariantMatrix::VariantContigMatrix {
                samples,
                target_names,
                target_lengths,
                ..
            } => {
                // Writes out coverage and variance values in metabat format
                // NOTE: This does not perform the metabat adjusted coverage formula
                // If you would like to use that then you are better off using CoverM
                // https://github.com/wwood/CoverM
                match table_type {
                    ReadType::Short => {
                        let file_name = format!("{}/rosella_coverages.tsv", &output,);

                        let file_path = Path::new(&file_name);

                        let mut file_open = match File::create(file_path) {
                            Ok(file) => file,
                            Err(e) => {
                                println!("Cannot create file {:?}", e);
                                std::process::exit(1)
                            }
                        };

                        // Write sample names in headers
                        write!(file_open, "{}", "contigName").unwrap();
                        write!(file_open, "\t{}", "contigLen").unwrap();
                        write!(file_open, "\t{}", "totalAvgDepth").unwrap();

                        for sample in samples.short.iter() {
                            write!(file_open, "\t{}.bam\t{}.bam-var", sample, sample).unwrap();
                        }

                        write!(file_open, "\n").unwrap();

                        for (tid, contig) in target_names.iter() {
                            match samples.short_coverages.get(tid) {
                                Some(cov) => {
                                    // Write the contig name
                                    write!(file_open, "{}", contig).unwrap();

                                    // Wirte the length
                                    let length = target_lengths.get(tid).unwrap();
                                    write!(file_open, "\t{}", length).unwrap();

                                    let var = samples.short_variances.get(tid).unwrap();

                                    // calc total avg depth from coverage values and write it to file
                                    let tot_avg_depth = cov.iter().sum::<f64>() / cov.len() as f64;
                                    write!(file_open, "\t{}", tot_avg_depth).unwrap();

                                    for (cov, var) in izip!(cov.iter(), var.iter()) {
                                        write!(file_open, "\t{}\t{}", cov, var).unwrap();
                                    }
                                    write!(file_open, "\n").unwrap();
                                }
                                None => {} // skip
                            }
                        }
                    }
                    ReadType::Long => {
                        let file_name = format!("{}/rosella_long_coverages.tsv", &output,);

                        let file_path = Path::new(&file_name);

                        let mut file_open = match File::create(file_path) {
                            Ok(file) => file,
                            Err(e) => {
                                println!("Cannot create file {:?}", e);
                                std::process::exit(1)
                            }
                        };

                        // Write sample names in headers
                        write!(file_open, "{}", "contigName").unwrap();
                        write!(file_open, "\t{}", "contigLen").unwrap();
                        write!(file_open, "\t{}", "totalAvgDepth").unwrap();

                        for sample in samples.long.iter() {
                            write!(file_open, "\t{}.bam\t{}.bam-var", sample, sample).unwrap();
                        }

                        write!(file_open, "\n").unwrap();

                        for (tid, contig) in target_names.iter() {
                            match samples.long_coverages.get(tid) {
                                Some(cov) => {
                                    // Write the contig name
                                    write!(file_open, "{}", contig).unwrap();

                                    // Wirte the length
                                    let length = target_lengths.get(tid).unwrap();
                                    write!(file_open, "\t{}", length).unwrap();

                                    let var = samples.long_variances.get(tid).unwrap();

                                    // calc total avg depth from coverage values and write it to file
                                    let tot_avg_depth = cov.iter().sum::<f64>() / cov.len() as f64;
                                    write!(file_open, "\t{}", tot_avg_depth).unwrap();

                                    for (cov, var) in izip!(cov.iter(), var.iter()) {
                                        write!(file_open, "\t{}\t{}", cov, var).unwrap();
                                    }
                                    write!(file_open, "\n").unwrap();
                                }
                                None => {} // skip
                            }
                        }
                    }
                    _ => panic!("ReadType {:?} not supported", table_type),
                }
            }
        }
    }

    fn merge_matrices(
        &mut self,
        tid: u32,
        other_matrix: VariantMatrix,
        assembly_sample_count: usize,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut variant_rates,
                ref mut samples,
                ref mut target_names,
                ref mut target_lengths,
                ref mut target_ids,
                ..
            } => {
                // Retrieve the needed parts from the other matrix
                let (other_rates, other_names, other_lengths, other_ids, other_samples) =
                    match other_matrix {
                        VariantMatrix::VariantContigMatrix {
                            variant_rates,
                            target_names,
                            target_lengths,
                            target_ids,
                            samples,
                            ..
                        } => (
                            variant_rates,
                            target_names,
                            target_lengths,
                            target_ids,
                            samples,
                        ),
                    };

                debug!("Other rates {:?}", &other_rates);
                // Trim out the assembly sample from each merge
                for (tid, rates) in other_rates.into_iter() {
                    let mut contig_rates = variant_rates.entry(tid).or_insert(HashMap::new());
                    for (var, var_rate) in rates.into_iter() {
                        contig_rates.insert(
                            var,
                            var_rate
                                .into_iter()
                                .take(assembly_sample_count)
                                .collect_vec(),
                        );
                    }
                }
                debug!("Other samples {:?}", &other_samples);

                for (tid, cov) in other_samples.short_coverages.iter() {
                    samples.short_coverages.insert(
                        *tid,
                        cov.iter()
                            .cloned()
                            .take(assembly_sample_count)
                            .collect_vec(),
                    );
                }

                for (tid, cov) in other_samples.long_coverages.iter() {
                    samples.long_coverages.insert(
                        *tid,
                        cov.iter()
                            .cloned()
                            .take(assembly_sample_count)
                            .collect_vec(),
                    );
                }

                for (tid, vari) in other_samples.short_variances.into_iter() {
                    samples.short_variances.insert(
                        tid,
                        vari.iter()
                            .cloned()
                            .take(assembly_sample_count)
                            .collect_vec(),
                    );
                }

                for (tid, vari) in other_samples.long_variances.into_iter() {
                    samples.long_variances.insert(
                        tid,
                        vari.iter()
                            .cloned()
                            .take(assembly_sample_count)
                            .collect_vec(),
                    );
                }

                debug!("Other names {:?}", &other_names);

                for (tid, name) in other_names.into_iter() {
                    target_names.insert(tid, name);
                }

                debug!("Other lengths {:?}", &other_lengths);

                for (tid, length) in other_lengths.into_iter() {
                    target_lengths.insert(tid, length);
                }

                debug!("Other lengths {:?}", &other_ids);

                for (name, tid) in other_ids.into_iter() {
                    target_ids.insert(name, tid);
                }

                if samples.short.contains(&"".to_string()) {
                    samples.short = other_samples
                        .short
                        .into_iter()
                        .take(assembly_sample_count)
                        .collect_vec();
                };

                if samples.long.contains(&"".to_string()) {
                    samples.long = other_samples
                        .long
                        .into_iter()
                        .take(assembly_sample_count)
                        .collect_vec();
                };
            }
        }
    }

    fn bin_contigs(&self, output: &str, m: &clap::ArgMatches) {
        match self {
            VariantMatrix::VariantContigMatrix { samples, .. } => {
                external_command_checker::check_for_flight();
                // let min_cluster_size
                let mut cmd_string;

                if samples.short.len() > 0 && samples.long.len() > 0 {
                    cmd_string = format!(
                        "flight bin \
                        --assembly {} \
                        --input {} \
                        --long_input {} \
                        --kmer_frequencies {} \
                        --min_contig_size {} \
                        --min_bin_size {} \
                        --n_neighbors {} \
                        --output_directory {} \
                        --cores {} ",
                        m.value_of("reference").unwrap(),
                        match m.is_present("coverage-values") {
                            true => m.value_of("coverage-values").unwrap().to_string(),
                            false => format!("{}/rosella_coverages.tsv", &output),
                        },
                        match m.is_present("longread-coverage-values") {
                            true => m.value_of("longread-coverage-values").unwrap().to_string(),
                            false => format!("{}/rosella_long_coverages.tsv", &output),
                        },
                        match m.is_present("kmer-frequencies") {
                            true => m.value_of("kmer-frequencies").unwrap().to_string(),
                            false => format!("{}/rosella_kmer_table.tsv", &output),
                        },
                        m.value_of("min-contig-size").unwrap(),
                        m.value_of("min-bin-size").unwrap(),
                        m.value_of("n-neighbors").unwrap(),
                        format!("{}/", &output),
                        m.value_of("threads").unwrap(),
                    );
                } else if samples.short.len() > 0 {
                    cmd_string = format!(
                        "flight bin \
                        --assembly {} \
                        --input {} \
                        --kmer_frequencies {} \
                        --min_contig_size {} \
                        --min_bin_size {} \
                        --n_neighbors {} \
                        --output_directory {} \
                        --cores {} ",
                        m.value_of("reference").unwrap(),
                        match m.is_present("coverage-values") {
                            true => m.value_of("coverage-values").unwrap().to_string(),
                            false => format!("{}/rosella_coverages.tsv", &output),
                        },
                        match m.is_present("kmer-frequencies") {
                            true => m.value_of("kmer-frequencies").unwrap().to_string(),
                            false => format!("{}/rosella_kmer_table.tsv", &output),
                        },
                        m.value_of("min-contig-size").unwrap(),
                        m.value_of("min-bin-size").unwrap(),
                        m.value_of("n-neighbors").unwrap(),
                        format!("{}/", &output),
                        m.value_of("threads").unwrap(),
                    );
                } else {
                    cmd_string = format!(
                        "flight bin \
                    --assembly {} \
                    --long_input {} \
                    --kmer_frequencies {} \
                    --min_contig_size {} \
                    --min_bin_size {} \
                    --n_neighbors {} \
                    --output_directory {} \
                    --cores {} ",
                        m.value_of("reference").unwrap(),
                        match m.is_present("longread-coverage-values") {
                            true => m.value_of("longread-coverage-values").unwrap().to_string(),
                            false => format!("{}/rosella_long_coverages.tsv", &output),
                        },
                        match m.is_present("kmer-frequencies") {
                            true => m.value_of("kmer-frequencies").unwrap().to_string(),
                            false => format!("{}/rosella_kmer_table.tsv", &output),
                        },
                        m.value_of("min-contig-size").unwrap(),
                        m.value_of("min-bin-size").unwrap(),
                        m.value_of("n-neighbors").unwrap(),
                        format!("{}/", &output),
                        m.value_of("threads").unwrap(),
                    );
                }

                command::finish_command_safely(
                    std::process::Command::new("bash")
                        .arg("-c")
                        .arg(&cmd_string)
                        .stderr(std::process::Stdio::null())
                        // .stdout(std::process::Stdio::piped())
                        .spawn()
                        .expect("Unable to execute bash"),
                    "flight",
                );
            }
        }
    }


    fn refine_bins(&self, output: &str, m: &clap::ArgMatches) {
        match self {
            VariantMatrix::VariantContigMatrix { samples, .. } => {
                external_command_checker::check_for_flight();
                // let min_cluster_size
                let mut cmd_string;

                if samples.short.len() > 0 && samples.long.len() > 0 {
                    cmd_string = format!(
                        "flight refine \
                            {} \
                            --extension {} \
                            --input {} \
                            --long_input {} \
                            --kmer_frequencies {} \
                            --min_contig_size {} \
                            --min_bin_size {} \
                            --n_neighbors {} \
                            --output_directory {} \
                            --cores {} --max_contamination {} {} {}",
                        if m.is_present("genome-fasta-directory") {
                            format!("--genome_directory {}", m.value_of("genome-fasta-directory").unwrap())
                        } else {
                            format!("--genome_paths {}", m.values_of("genome-fasta-files").unwrap().collect::<Vec<&str>>().join(" "))
                        },
                        m.value_of("genome-fasta-extension").unwrap(),
                        match m.is_present("coverage-values") {
                            true => m.value_of("coverage-values").unwrap().to_string(),
                            false => format!("{}/rosella_coverages.tsv", &output),
                        },
                        match m.is_present("longread-coverage-values") {
                            true => m.value_of("longread-coverage-values").unwrap().to_string(),
                            false => format!("{}/rosella_long_coverages.tsv", &output),
                        },
                        match m.is_present("kmer-frequencies") {
                            true => m.value_of("kmer-frequencies").unwrap().to_string(),
                            false => format!("{}/rosella_kmer_table.tsv", &output),
                        },
                        m.value_of("min-contig-size").unwrap(),
                        m.value_of("min-bin-size").unwrap(),
                        m.value_of("n-neighbors").unwrap(),
                        format!("{}/", &output),
                        m.value_of("threads").unwrap(),
                        m.value_of("max-contamination").unwrap(),
                        if m.is_present("checkm-file") {
                            format!("--checkm_file {}", m.value_of("checkm-file").unwrap())
                        } else {
                            "".to_string()
                        },
                        if m.is_present("contaminated-only") {
                            "--contaminated_only True"
                        } else {
                            ""
                        },
                    );
                } else if samples.short.len() > 0 {
                    cmd_string = format!(
                        "flight refine \
                            {} \
                            --extension {} \
                            --input {} \
                            --kmer_frequencies {} \
                            --min_contig_size {} \
                            --min_bin_size {} \
                            --n_neighbors {} \
                            --output_directory {} \
                            --cores {} --max_contamination {} {} {}",
                        if m.is_present("genome-fasta-directory") {
                            format!("--genome_directory {}", m.value_of("genome-fasta-directory").unwrap())
                        } else {
                            format!("--genome_paths {}", m.values_of("genome-fasta-files").unwrap().collect::<Vec<&str>>().join(" "))
                        },
                        m.value_of("genome-fasta-extension").unwrap(),
                        match m.is_present("coverage-values") {
                            true => m.value_of("coverage-values").unwrap().to_string(),
                            false => format!("{}/rosella_coverages.tsv", &output),
                        },
                        match m.is_present("kmer-frequencies") {
                            true => m.value_of("kmer-frequencies").unwrap().to_string(),
                            false => format!("{}/rosella_kmer_table.tsv", &output),
                        },
                        m.value_of("min-contig-size").unwrap(),
                        m.value_of("min-bin-size").unwrap(),
                        m.value_of("n-neighbors").unwrap(),
                        format!("{}/", &output),
                        m.value_of("threads").unwrap(),
                        m.value_of("max-contamination").unwrap(),
                        if m.is_present("checkm-file") {
                            format!("--checkm_file {}", m.value_of("checkm-file").unwrap())
                        } else {
                            "".to_string()
                        },
                        if m.is_present("contaminated-only") {
                            "--contaminated_only True"
                        } else {
                            ""
                        },
                    );
                } else {
                    cmd_string = format!(
                        "flight refine \
                            {} \
                            --extension {} \
                            --long_input {} \
                            --kmer_frequencies {} \
                            --min_contig_size {} \
                            --min_bin_size {} \
                            --n_neighbors {} \
                            --output_directory {} \
                            --cores {} --max_contamination {} {} {}",
                        if m.is_present("genome-fasta-directory") {
                            format!("--genome_directory {}", m.value_of("genome-fasta-directory").unwrap())
                        } else {
                            format!("--genome_paths {}", m.values_of("genome-fasta-files").unwrap().collect::<Vec<&str>>().join(" "))
                        },
                        m.value_of("genome-fasta-extension").unwrap(),
                        match m.is_present("longread-coverage-values") {
                            true => m.value_of("longread-coverage-values").unwrap().to_string(),
                            false => format!("{}/rosella_long_coverages.tsv", &output),
                        },
                        match m.is_present("kmer-frequencies") {
                            true => m.value_of("kmer-frequencies").unwrap().to_string(),
                            false => format!("{}/rosella_kmer_table.tsv", &output),
                        },
                        m.value_of("min-contig-size").unwrap(),
                        m.value_of("min-bin-size").unwrap(),
                        m.value_of("n-neighbors").unwrap(),
                        format!("{}/", &output),
                        m.value_of("threads").unwrap(),
                        m.value_of("max-contamination").unwrap(),
                        if m.is_present("checkm-file") {
                            format!("--checkm_file {}", m.value_of("checkm-file").unwrap())
                        } else {
                            "".to_string()
                        },
                        if m.is_present("contaminated-only") {
                            "--contaminated_only True"
                        } else {
                            ""
                        },
                    );
                }

                command::finish_command_safely(
                    std::process::Command::new("bash")
                        .arg("-c")
                        .arg(&cmd_string)
                        .stderr(std::process::Stdio::null())
                        // .stdout(std::process::Stdio::piped())
                        .spawn()
                        .expect("Unable to execute bash"),
                    "flight",
                );
            }
        }
    }

    fn finalize_bins(&mut self, output: &str, m: &clap::ArgMatches) {
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                target_lengths,
                samples,
                ref mut final_bins,
                ..
            } => {
                let min_bin_size = m.value_of("min-bin-size").unwrap().parse().unwrap();
                let min_contig_size = m.value_of("min-contig-size").unwrap().parse().unwrap();

                let mut bins =
                    std::fs::File::open(format!("{}/rosella_bins.json", &output,)).unwrap();
                let mut data = String::new();
                bins.read_to_string(&mut data).unwrap();
                let mut bins: HashMap<usize, Vec<i32>> = serde_json::from_str(&data).unwrap();

                // Ordered vectors for large and small contigs. Should be identical to python
                // let mut large_contigs = Vec::new();
                // let mut small_contigs = Vec::new();

                let small_contigs: Vec<i32> = target_names
                    .par_iter()
                    .filter_map(|(tid, _)| {
                        if &1000 <= target_lengths.get(tid).unwrap()
                            && target_lengths.get(tid).unwrap() < &min_contig_size
                        {
                            Some(*tid)
                        } else {
                            None
                        }
                    })
                    .collect();

                // Remove small bins and collect contigs
                let mut removed_contigs: Vec<i32> = Vec::new();
                let mut removed_bins: Vec<usize> = Vec::new();
                for (bin, contigs) in bins.iter() {
                    if bin == &0 {
                        removed_bins.push(*bin);
                        removed_contigs.par_extend(contigs.par_iter());
                    } else {
                        let bin_length: u64 = contigs
                            .par_iter()
                            .fold_with(0, |acc, tid| {
                                // let tid = large_contigs[id];
                                let length = target_lengths.get(&tid).unwrap();
                                acc + *length as u64
                            })
                            .sum::<u64>();
                        if bin_length < min_bin_size {
                            // Break up the bin and send off the contigs to other bins
                            removed_bins.push(*bin);
                            removed_contigs.par_extend(contigs.par_iter());
                            // contigs.par_iter().for_each(|idx| {
                            //
                            // });
                        }
                    }
                }

                // remove bins
                for bin in removed_bins {
                    bins.remove(&bin).unwrap();
                }

                // If we have enough samples, we can attempt to rescue the contigs we just
                // removed.

                if samples.short.len() >= 3 || samples.long.len() >= 3 {
                    // Iterate through removed contigs and check their average correlation with
                    // other bins
                    // correlate_with_bins(&removed_contigs[..], &mut bins, &samples, &target_lengths);

                    // Now we use a similar method to above to rescue small contigs
                    // correlate_with_bins(&small_contigs[..], &mut bins, &samples, &target_lengths);
                }

                *final_bins = bins
            }
        }
    }

    fn write_bins(&self, output: &str, reference: &str, mode: &str) {
        match self {
            VariantMatrix::VariantContigMatrix {
                final_bins,
                target_names,
                ..
            } => {
                let mut reference_file = generate_faidx(reference);
                final_bins.iter().for_each(|(bin, contigs)| {

                    let mut writer = bio::io::fasta::Writer::to_file(format!(
                        "{}/rosella_{}.{}.fna",
                        &output, mode, bin
                    ))
                    .unwrap();
                    let mut ref_seq = Vec::new();

                    for tid in contigs.iter() {
                        // Update all contig information
                        fetch_contig_from_reference(
                            &mut reference_file,
                            &target_names.get(tid).unwrap().as_bytes().to_vec(),
                        );
                        read_sequence_to_vec(
                            &mut ref_seq,
                            &mut reference_file,
                            &target_names.get(tid).unwrap().as_bytes().to_vec(),
                        );
                        writer
                            .write(&target_names.get(tid).unwrap(), None, &ref_seq[..])
                            .expect("Unable to write FASTA record");
                    }
                });
            }
        }
    }

    fn read_inputs(
        &mut self,
        coverage_path: Option<&str>,
        long_coverage_path: Option<&str>,
        kmer_path: Option<&str>,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                target_lengths,
                target_ids,
                ref mut samples,
                ref mut variant_rates,
                ref mut kfrequencies,
                ..
            } => {
                match coverage_path {
                    Some(coverage_str) => {
                        // Read the coverage information in
                        let coverage_file = std::fs::File::open(coverage_str).unwrap();
                        let coverage_buffer = std::io::BufReader::new(coverage_file);
                        for (line_idx, line) in coverage_buffer.lines().enumerate() {
                            let line = line.unwrap();
                            let values = line.split('\t').collect_vec();

                            if line_idx == 0 {
                                for value in values[3..].iter().step_by(2) {
                                    let value = value.replace(".bam", "");
                                    samples.short.push(value);
                                }
                            } else {
                                let tid = match target_ids.get(values[0]) {
                                    Some(v) => v,
                                    None => panic!("Could not find contig: {:?} in assembly", values[0])
                                };

                                for (idx, value) in values[3..].iter().enumerate() {
                                    if idx % 2 == 0 {
                                        let idx = idx / 2;
                                        let contig_coverages = samples
                                            .short_coverages
                                            .entry(*tid)
                                            .or_insert(vec![0.; samples.short.len()]);
                                        contig_coverages[idx as usize] = value.parse().unwrap();
                                    } else {
                                        let idx = idx / 2; // Rust defaults to floor division, so this should yield correct index
                                        let contig_variances = samples
                                            .short_variances
                                            .entry(*tid)
                                            .or_insert(vec![0.; samples.short.len()]);
                                        contig_variances[idx as usize] = value.parse().unwrap();
                                    }
                                }
                            }
                        }
                    }
                    None => {
                        // Nothing to read, move on
                    }
                }

                match long_coverage_path {
                    Some(coverage_str) => {
                        // Read the coverage information in
                        // Read the coverage information in
                        let coverage_file = std::fs::File::open(coverage_str).unwrap();

                        let coverage_buffer = std::io::BufReader::new(coverage_file);
                        for (line_idx, line) in coverage_buffer.lines().enumerate() {
                            let line = line.unwrap();
                            let values = line.split('\t').collect_vec();

                            if line_idx == 0 {
                                for value in values[3..].iter().step_by(2) {
                                    let value = value.replace(".bam", "");
                                    samples.long.push(value);
                                }
                            } else {
                                let tid = match target_ids.get(values[0]) {
                                    Some(v) => v,
                                    None => panic!("Could not find contig: {:?} in assembly", values[0])
                                };

                                for (idx, value) in values[3..].iter().enumerate() {
                                    if idx % 2 == 0 {
                                        let idx = idx / 2;
                                        let contig_coverages = samples
                                            .long_coverages
                                            .entry(*tid)
                                            .or_insert(vec![0.; samples.long.len()]);
                                        contig_coverages[idx as usize] = value.parse().unwrap();
                                    } else {
                                        let idx = idx / 2; // Rust defaults to floor division, so this should yield correct index
                                        let contig_variances = samples
                                            .long_variances
                                            .entry(*tid)
                                            .or_insert(vec![0.; samples.long.len()]);
                                        contig_variances[idx as usize] = value.parse().unwrap();
                                    }
                                }
                            }
                        }
                    }
                    None => {
                        // Nothing to read, move on
                    }
                }

                match kmer_path {
                    Some(kmer_str) => {
                        // Read the kmer information in
                        let kmer_file = std::fs::File::open(kmer_str).unwrap();

                        let kmer_buffer = std::io::BufReader::new(kmer_file);
                        for (line_idx, line) in kmer_buffer.lines().enumerate() {
                            let line = line.unwrap();
                            let values = line.split('\t').collect_vec();

                            if line_idx == 0 {
                                for value in values[2..].iter() {
                                    kfrequencies
                                        .entry(value.as_bytes().to_vec())
                                        .or_insert(vec![0; target_names.len()]);
                                }
                            } else {
                                let tid = match target_ids.get(values[0]) {
                                    Some(v) => v,
                                    None => panic!("Could not find contig: {:?} in assembly", values[0])
                                };

                                for (value, (_kmer, kmer_vec)) in
                                    values[2..].iter().zip(kfrequencies.iter_mut())
                                {
                                    kmer_vec[*tid as usize] = value.parse().unwrap();
                                }
                            }
                        }
                    }
                    None => {
                        // Nothing to read, move on
                    }
                }
            }
        }
    }
}

/// Add read count entry to cluster hashmap
pub fn add_entry(
    shared_read_counts: &mut HashMap<usize, HashMap<usize, usize>>,
    clust1: usize,
    clust2: usize,
    count: usize,
) {
    let entry = shared_read_counts.entry(clust1).or_insert(HashMap::new());
    entry.insert(clust2, count);
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::variants;
    use std::collections::HashSet;
}
