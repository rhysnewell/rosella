use crate::external_command_checker;
use bio::alignment::sparse::hash_kmers;
use bio::io::fasta::IndexedReader;
use bio::stats::{
    bayesian,
    probs::{LogProb, PHREDProb, Prob},
};
use bird_tool_utils::command;
use estimation::contig_variants::*;
use itertools::izip;
use itertools::Itertools;
use model::variants::*;
use ordered_float::NotNan;
use rayon::prelude::*;
use rgsl::statistics::spearman;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::str;
use std::sync::mpsc::channel;
use utils::{fetch_contig_from_reference, generate_faidx, read_sequence_to_vec};

#[derive(Debug, Clone)]
/// Container for all variants within a genome and associated clusters
pub enum VariantMatrix<'b> {
    VariantContigMatrix {
        coverages: HashMap<i32, Vec<f64>>,
        average_genotypes: HashMap<i32, Vec<f64>>,
        variances: HashMap<i32, Vec<f64>>,
        // TID, Position, Variant Type, Base
        all_variants: HashMap<i32, HashMap<i64, HashMap<Variant, Base>>>,
        // Placeholder hashmap for the depths of each contig for a sample
        // Deleted after use
        depths: HashMap<i32, Vec<i32>>,
        target_names: BTreeMap<i32, String>,
        target_lengths: HashMap<i32, f64>,
        sample_names: Vec<String>,
        kfrequencies: BTreeMap<Vec<u8>, Vec<usize>>,
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
    pub fn new_matrix(sample_count: usize) -> VariantMatrix<'b> {
        VariantMatrix::VariantContigMatrix {
            coverages: HashMap::new(),
            variances: HashMap::new(),
            average_genotypes: HashMap::new(),
            all_variants: HashMap::new(),
            depths: HashMap::new(),
            target_names: BTreeMap::new(),
            target_lengths: HashMap::new(),
            sample_names: vec!["".to_string(); sample_count],
            kfrequencies: BTreeMap::new(),
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
    fn get_sample_name(&self, sample_idx: usize) -> &str;

    /// Returns the total amount of variant alleles
    fn get_variant_count(&self) -> i64;

    fn add_sample_name(&mut self, sample_name: String, sample_idx: usize);

    fn add_info(&mut self, tid: usize, target_name: Vec<u8>, target_len: u64);

    /// Adds the variant information retrieved from VCF records on a per reference basis
    fn add_variant_to_matrix(&mut self, sample_idx: usize, base: &Base, tid: usize);

    /// Remove variants that were detected outside of suitable coverage values
    fn remove_variants(&mut self, tid: i32, sample_idx: usize, contig_info: Vec<f64>);

    /// Per sample FDR calculation
    fn remove_false_discoveries(&mut self, tid: i32, alpha: f64, reference: &str);

    /// Returns the alleles at the current position
    /// as a mutable reference
    fn variants(&mut self, tid: i32, pos: i64) -> Option<&mut HashMap<Variant, Base>>;

    /// Returns all variants found within a contig
    fn variants_of_contig(&mut self, tid: i32)
        -> Option<&mut HashMap<i64, HashMap<Variant, Base>>>;

    /// Takes [VariantStats](contig_variants/VariantStats) struct for single contig and adds to
    /// [VariantMatrix](VariantMatrix)
    fn add_contig(&mut self, variant_stats: VariantStats, sample_count: usize, sample_idx: usize);

    /// Calculates the kmer frequencies for a given contig
    fn calc_kmer_frequencies(
        &mut self,
        tid: u32,
        kmer_size: usize,
        reference_file: &mut IndexedReader<File>,
        contig_count: usize,
    );

    fn write_kmer_table(&self, output: &str);

    fn calc_variant_rates(&mut self, tid: u32, window_size: usize, sample_idx: usize);

    fn write_variant_rates(&self, output: &str);

    fn write_coverage(&self, output: &str);

    fn merge_matrices(
        &mut self,
        tid: u32,
        other_matrix: VariantMatrix,
        assembly_sample_count: usize,
    );

    fn bin_contigs(&self, output: &str, m: &clap::ArgMatches);

    fn finalize_bins(&mut self, output: &str, m: &clap::ArgMatches);

    fn write_bins(&self, output: &str, reference: &str);

    fn read_inputs(
        &mut self,
        coverages_path: Option<&str>,
        kmer_path: Option<&str>,
        rates_path: Option<&str>,
    );
}

#[allow(unused)]
impl VariantMatrixFunctions for VariantMatrix<'_> {
    fn setup(&mut self) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut coverages,
                ref mut average_genotypes,
                ref mut all_variants,
                ref mut target_names,
                ref mut target_lengths,
                ref mut sample_names,
                ref mut kfrequencies,
                ..
            } => {
                *coverages = HashMap::new();
                *average_genotypes = HashMap::new();
                *all_variants = HashMap::new();
                *target_names = BTreeMap::new();
                *target_lengths = HashMap::new();
                *sample_names = vec![];
                *kfrequencies = BTreeMap::new();
            }
        }
    }

    fn get_variant_count(&self) -> i64 {
        match self {
            VariantMatrix::VariantContigMatrix { variant_counts, .. } => {
                let mut total = 0;
                for (_, count) in variant_counts {
                    total += count;
                }
                total as i64
            }
        }
    }

    fn get_sample_name(&self, sample_idx: usize) -> &str {
        match self {
            VariantMatrix::VariantContigMatrix { sample_names, .. } => {
                if sample_names[sample_idx].to_string().contains(".tmp") {
                    return &sample_names[sample_idx][15..];
                } else {
                    return &sample_names[sample_idx];
                }
            }
        }
    }

    fn add_sample_name(&mut self, sample_name: String, sample_idx: usize) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut sample_names,
                ..
            } => {
                debug!("adding sample {} at index {}", &sample_name, &sample_idx);
                sample_names[sample_idx] = sample_name;
            }
        }
    }

    fn add_info(&mut self, tid: usize, target_name: Vec<u8>, target_len: u64) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut sample_names,
                ref mut target_names,
                ref mut target_lengths,
                ..
            } => {
                target_names
                    .entry(tid as i32)
                    .or_insert(String::from_utf8(target_name).unwrap()); // add contig name

                target_lengths
                    .entry(tid as i32)
                    .or_insert(target_len as f64); // add contig length
            }
        }
    }

    fn add_variant_to_matrix(&mut self, sample_idx: usize, base: &Base, tid: usize) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut sample_names,
                ref mut all_variants,
                ref mut target_names,
                ref mut target_lengths,
                ref mut variant_counts,
                ..
            } => {
                // Initialize contig id in variant hashmap
                let contig_variants = all_variants.entry(tid as i32).or_insert(HashMap::new());

                let mut con_var_counts = variant_counts.entry(tid as i32).or_insert(0);
                // Apppend the sample index to each variant abundance
                // Initialize the variant position index
                // Also turns out to be the total number of variant positions
                let position_variants = contig_variants.entry(base.pos).or_insert(HashMap::new());
                let allele = position_variants
                    .entry(base.variant.clone())
                    .or_insert(base.clone());

                allele.combine_sample(base, sample_idx, 0);
                *con_var_counts += 1; // increment count for contig
            }
        }
    }

    fn remove_variants(&mut self, tid: i32, sample_idx: usize, stats: Vec<f64>) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut all_variants,
                ..
            } => {
                let upper_limit = stats[0] + stats[1];
                let lower_limit = stats[0] - stats[1];
                match all_variants.get_mut(&tid) {
                    Some(contig_variants) => {
                        let (to_remove_s, to_remove_r) = std::sync::mpsc::channel();

                        contig_variants.par_iter().for_each_with(
                            to_remove_s,
                            |s, (position, variants)| {
                                let mut total_depth = 0;
                                for (_, base) in variants {
                                    total_depth += base.depth[sample_idx];
                                }
                                // Scan within 10bp of variant to see if there is
                                // >= 3 SNVs
                                let lower_pos = std::cmp::min(position - 10, 0);
                                let upper_pos = position + 10;
                                let (count_s, count_r) = std::sync::mpsc::channel();
                                (lower_pos..upper_pos).into_par_iter().for_each_with(
                                    count_s,
                                    |s_2, window| {
                                        if &window != position {
                                            if contig_variants.contains_key(&window) {
                                                let window_variants =
                                                    contig_variants.get(&window).unwrap();
                                                for (window_variant, _) in window_variants.iter() {
                                                    match window_variant {
                                                        &Variant::SNV(_) => s_2.send(1),
                                                        _ => s_2.send(0),
                                                    };
                                                }
                                            }
                                        }
                                    },
                                );
                                let count: i64 = count_r.iter().sum();
                                let count = 0;
                                let total_depth = total_depth as f64;
                                if (total_depth < lower_limit || total_depth > upper_limit)
                                    && total_depth > 0.
                                    || count >= 3
                                {
                                    s.send(*position);
                                }
                            },
                        );
                        let to_remove: Vec<i64> = to_remove_r.iter().collect();

                        for position in to_remove.iter() {
                            contig_variants.remove_entry(position).unwrap();
                        }
                    }
                    None => {
                        // do nothing
                    }
                }
            }
        }
    }

    fn remove_false_discoveries(&mut self, tid: i32, alpha: f64, reference: &str) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut sample_names,
                ref mut all_variants,
                ref mut target_names,
                ref mut target_lengths,
                ..
            } => {
                // collect variant probabilities e.g. PHRED Sums
                let mut prob_dist = Vec::new();
                for (tid, position_map) in all_variants.iter() {
                    position_map.iter().for_each(|(pos, alleles)| {
                        for (variant, allele) in alleles.iter() {
                            if allele.variant != Variant::None {
                                // get the mean PHRED value
                                // This is probably not the best method
                                let mut allele_probs = allele.quals.par_iter().sum();
                                // allele_probs.par_sort();

                                // let mut allele_probs = allele_probs
                                //     .iter()
                                //     .map(|p| {
                                //         LogProb::from(Prob(10_f64.powf(-(**p as f64) / 10)))
                                //     })
                                //     .collect_vec();

                                prob_dist.push(
                                    NotNan::new(allele_probs).expect("Unable to convert to NotNan"),
                                );
                            }
                        }
                    });
                }
                // sort prob_dist
                prob_dist.par_sort();
                // turn into logprob values
                let prob_dist = prob_dist
                    .into_iter()
                    .rev()
                    .map(|p| LogProb::from(PHREDProb(*p)))
                    .collect_vec();

                // info!("Prob {:?}", &prob_dist);

                // turn prob dist into posterior exact probabilities
                debug!("Calculating PEP dist...");
                let pep_dist = prob_dist
                    .iter()
                    .map(|p| p.ln_one_minus_exp())
                    .collect::<Vec<LogProb>>();

                // info!("PEP {:?}", &pep_dist);

                // calculate fdr values
                let fdrs = bayesian::expected_fdr(&pep_dist);
                // info!("FDRs {:?}", &fdrs);

                let alpha = LogProb::from(Prob(alpha));
                let mut threshold = None;
                if fdrs.is_empty() {
                    debug!("Genome {}: FDR calculations are empty...", reference);
                    threshold = None;
                } else if fdrs[0] > alpha {
                    debug!(
                        "Genome {}: FDR threshold for alpha of {:?} calculated as {:?}",
                        reference,
                        alpha,
                        LogProb::ln_one()
                    );
                    threshold = Some(LogProb::ln_one());
                } else {
                    // find the largest pep for which fdr <= alpha
                    // do not let peps with the same value cross the boundary
                    for i in (0..fdrs.len()).rev() {
                        if fdrs[i] <= alpha && (i == 0 || pep_dist[i] != pep_dist[i - 1]) {
                            let prob = prob_dist[i];
                            debug!(
                                "Genome {}: FDR threshold for alpha of {:?} calculated as {:?}",
                                reference, alpha, &prob
                            );
                            threshold = Some(prob);
                            break;
                        }
                    }
                }
                match threshold {
                    Some(prob) => {
                        let mut total_removed = 0;
                        let mut total_kept = 0;
                        for (tid, position_map) in all_variants.iter_mut() {
                            let mut positions_to_remove = Vec::new();
                            for (pos, alleles) in position_map.iter_mut() {
                                let mut to_remove = Vec::new();
                                let mut variant_alleles = 0;
                                for (variant, allele) in alleles.iter_mut() {
                                    if allele.variant != Variant::None {
                                        let mut allele_probs = allele.quals.par_iter().sum();
                                        let sum_prob = LogProb::from(PHREDProb(allele_probs));
                                        if sum_prob > prob {
                                            to_remove.push(variant.clone());
                                            total_removed += 1;
                                        } else {
                                            total_kept += 1;
                                        }
                                        variant_alleles += 1;
                                    }
                                }
                                if to_remove.len() > 1 {
                                    if to_remove.len() == variant_alleles {
                                        positions_to_remove.push(*pos);
                                    } else {
                                        for removing in to_remove.iter() {
                                            alleles
                                                .remove_entry(removing)
                                                .expect("Unable to remove variant entry");
                                        }
                                    }
                                }
                            }
                            if positions_to_remove.len() > 0 {
                                for position in positions_to_remove.iter() {
                                    position_map
                                        .remove_entry(position)
                                        .expect("Unable to remove position");
                                }
                            }
                        }
                        debug!("Genome {}: {} variants passed FDR threshold, {} did not pass FDR threshold", reference, total_kept, total_removed);
                        *all_variants = all_variants.clone();
                    }
                    None => {}
                }
            }
        }
    }

    fn variants(&mut self, tid: i32, pos: i64) -> Option<&mut HashMap<Variant, Base>> {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut all_variants,
                ..
            } => match all_variants.get_mut(&tid) {
                Some(contig_variants) => contig_variants.get_mut(&pos),
                None => None,
            },
        }
    }

    fn variants_of_contig(
        &mut self,
        tid: i32,
    ) -> Option<&mut HashMap<i64, HashMap<Variant, Base>>> {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut all_variants,
                ..
            } => all_variants.get_mut(&tid),
        }
    }

    fn add_contig(&mut self, variant_stats: VariantStats, sample_count: usize, sample_idx: usize) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut coverages,
                ref mut all_variants,
                target_names,
                ref mut variances,
                ref mut depths,
                ..
            } => match variant_stats {
                VariantStats::VariantContigStats {
                    tid,
                    coverage,
                    variance,
                    depth,
                    ..
                } => {
                    let var = variances
                        .entry(tid)
                        .or_insert(vec![0.0 as f64; sample_count]);
                    var[sample_idx] = variance;
                    let cov = coverages
                        .entry(tid)
                        .or_insert(vec![0.0 as f64; sample_count]);
                    cov[sample_idx] = coverage;

                    let contig_variants = all_variants.entry(tid).or_insert(HashMap::new());
                    for (pos, d) in depth.iter().enumerate() {
                        let position_variants =
                            contig_variants.entry(pos as i64).or_insert(HashMap::new());
                        let ref_depth = match position_variants.get(&Variant::None) {
                            Some(base) => base.truedepth[sample_idx],
                            None => 0,
                        };
                        for (variant, base_info) in position_variants.iter_mut() {
                            base_info.add_depth(sample_idx, *d, ref_depth);
                        }
                    }
                    depths.insert(tid, depth);
                }
            },
        }
    }

    fn calc_kmer_frequencies(
        &mut self,
        tid: u32,
        kmer_size: usize,
        reference_file: &mut IndexedReader<File>,
        contig_count: usize,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut kfrequencies,
                sample_names,
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
                let kmers = hash_kmers(&ref_seq[..], kmer_size);
                // Get kmer counts in a contig
                for (kmer, pos) in kmers {
                    let k = kfrequencies
                        .entry(kmer.to_vec())
                        .or_insert(vec![0; contig_count]);
                    // Insert kmer count at contig position
                    k[tid as usize] = pos.len();
                }
            }
        }
    }

    fn write_kmer_table(&self, output: &str) {
        match self {
            VariantMatrix::VariantContigMatrix {
                kfrequencies,
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
                for (kmer, _) in kfrequencies.iter() {
                    write!(file_open, "\t{}", str::from_utf8(&kmer[..]).unwrap()).unwrap();
                }

                write!(file_open, "\n").unwrap();

                for (idx, contig) in target_names.iter() {
                    write!(file_open, "{}", contig).unwrap();

                    // Write the length
                    let length = target_lengths.get(idx).unwrap();
                    write!(file_open, "\t{}", length).unwrap();

                    for (_kmer, counts) in kfrequencies.iter() {
                        write!(file_open, "\t{}", counts[*idx as usize]).unwrap();
                    }
                    write!(file_open, "\n").unwrap();
                }
            }
        }
    }

    fn calc_variant_rates(&mut self, tid: u32, window_size: usize, sample_idx: usize) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut all_variants,
                sample_names,
                target_names,
                target_lengths,
                ref mut variant_rates,
                ..
            } => {
                let placeholder_map = HashMap::new();
                let placeholder_variants = HashMap::new();
                let placeholder_name = "NOT_FOUND".to_string();
                let contig_variants = match all_variants.get(&(tid as i32)) {
                    Some(variants) => variants,
                    None => &placeholder_map,
                };
                let target_len = match target_lengths.get(&(tid as i32)) {
                    Some(length) => *length as usize,
                    None => 0,
                };
                let target_name = match target_names.get(&(tid as i32)) {
                    Some(name) => name,
                    None => &placeholder_name,
                };

                // If the contig size is small than the sliding window, then the variant rate
                // is knocked down to 0. This shouldn't matter so much since small contigs
                // are filtered out later on
                if target_len >= window_size {
                    let (rate_s, rate_r) = channel();

                    // Setup N sliding windows
                    (0..target_len - window_size)
                        .into_par_iter()
                        .for_each_with(rate_s, |s, i| {
                            // Collect all variants in current window
                            let mut counts = vec![0, 0];

                            (i..i + window_size).into_iter().for_each(|w_i| {
                                let pos_variants = match contig_variants.get(&(w_i as i64)) {
                                    Some(variants) => variants,
                                    None => &placeholder_variants,
                                };
                                for (variant, base) in pos_variants.iter() {
                                    match variant {
                                        Variant::SNV(_) => {
                                            if base.truedepth[sample_idx] > 0 {
                                                counts[0] += 1
                                            }
                                        }
                                        Variant::None => {}
                                        _ => {
                                            if base.truedepth[sample_idx] > 0 {
                                                counts[1] += 1
                                            }
                                        }
                                    }
                                }
                            });
                            // Calculate means for current window
                            let mut means = counts
                                .iter()
                                .map(|x| *x as f32 / window_size as f32)
                                .collect();
                            s.send(means).unwrap();
                        });

                    let n_windows = (target_len - window_size - 1) as f32;
                    // Sum the mean values of window variant rates, 0 is SNVs, 1 is SVs

                    let rates: Vec<Vec<f32>> = rate_r.iter().collect();
                    debug!("Rates {:?}", &rates);

                    let mut rates_sum = vec![0., 0.];
                    for rate in rates.iter() {
                        rates_sum[0] += rate[0];
                        rates_sum[1] += rate[1];
                    }

                    let rates: Vec<f32> = rates_sum
                        .iter()
                        .map(|rate| *rate / n_windows as f32)
                        .collect();

                    // Append rates to variant_rates variables
                    let mut contig_rates =
                        variant_rates.entry(tid as i32).or_insert(HashMap::new());
                    contig_rates
                        .entry(Var::SNV)
                        .or_insert(vec![0.; sample_names.len()])[sample_idx] = rates[0];
                    contig_rates
                        .entry(Var::SV)
                        .or_insert(vec![0.; sample_names.len()])[sample_idx] = rates[1];
                } else {
                    // set variant rates to zero for small contigs
                    let mut contig_rates =
                        variant_rates.entry(tid as i32).or_insert(HashMap::new());
                    contig_rates
                        .entry(Var::SNV)
                        .or_insert(vec![0.; sample_names.len()]);
                    contig_rates
                        .entry(Var::SV)
                        .or_insert(vec![0.; sample_names.len()]);
                }
            }
        }
    }

    fn write_variant_rates(&self, output: &str) {
        match self {
            VariantMatrix::VariantContigMatrix {
                sample_names,
                target_names,
                target_lengths,
                variant_rates,
                ..
            } => {
                let file_name = format!("{}/rosella_variant_rates.tsv", &output,);

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
                for sample in sample_names.iter() {
                    write!(file_open, "\t{}-SNV\t{}-SV", sample, sample).unwrap();
                }

                write!(file_open, "\n").unwrap();

                for (tid, contig) in target_names.iter() {
                    write!(file_open, "{}", contig).unwrap();

                    // Wirte the length
                    let length = target_lengths.get(tid).unwrap();
                    write!(file_open, "\t{}", length).unwrap();

                    let rates = variant_rates.get(tid).unwrap();
                    let snv_rates = rates.get(&Var::SNV).unwrap();
                    let sv_rates = rates.get(&Var::SV).unwrap();
                    for (snv, sv) in izip!(snv_rates.iter(), sv_rates.iter()) {
                        write!(file_open, "\t{}\t{}", snv, sv).unwrap();
                    }
                    write!(file_open, "\n").unwrap();
                }
            }
        }
    }

    fn write_coverage(&self, output: &str) {
        match self {
            VariantMatrix::VariantContigMatrix {
                sample_names,
                target_names,
                target_lengths,
                coverages,
                variances,
                ..
            } => {
                // Writes out coverage and variance values in metabat format
                // NOTE: This does not perform the metabat adjusted coverage formula
                // If you would like to use that then you are better off using CoverM
                // https://github.com/wwood/CoverM
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

                for sample in sample_names.iter() {
                    write!(file_open, "\t{}.bam\t{}.bam-var", sample, sample).unwrap();
                }

                write!(file_open, "\n").unwrap();

                for (tid, contig) in target_names.iter() {
                    // Write the contig name
                    write!(file_open, "{}", contig).unwrap();

                    // Wirte the length
                    let length = target_lengths.get(tid).unwrap();
                    write!(file_open, "\t{}", length).unwrap();

                    let cov = coverages.get(tid).unwrap();
                    let var = variances.get(tid).unwrap();

                    // calc total avg depth from coverage values and write it to file
                    let tot_avg_depth = cov.iter().sum::<f64>() / cov.len() as f64;
                    write!(file_open, "\t{}", tot_avg_depth).unwrap();

                    for (cov, var) in izip!(cov.iter(), var.iter()) {
                        write!(file_open, "\t{}\t{}", cov, var).unwrap();
                    }
                    write!(file_open, "\n").unwrap();
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
                ref mut coverages,
                ref mut variances,
                ref mut target_names,
                ref mut target_lengths,
                ref mut sample_names,
                ..
            } => {
                // Retrieve the needed parts from the other matrix
                let (other_rates, other_cov, other_var, other_names, other_lengths, other_samples) =
                    match other_matrix {
                        VariantMatrix::VariantContigMatrix {
                            variant_rates,
                            coverages,
                            variances,
                            target_names,
                            target_lengths,
                            sample_names,
                            ..
                        } => (
                            variant_rates,
                            coverages,
                            variances,
                            target_names,
                            target_lengths,
                            sample_names,
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

                for (tid, cov) in other_cov.into_iter() {
                    coverages.insert(
                        tid,
                        cov.into_iter().take(assembly_sample_count).collect_vec(),
                    );
                }

                for (tid, vari) in other_var.into_iter() {
                    variances.insert(
                        tid,
                        vari.into_iter().take(assembly_sample_count).collect_vec(),
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

                debug!("Other samples {:?}", &other_samples);
                if sample_names.contains(&"".to_string()) {
                    *sample_names = other_samples
                        .into_iter()
                        .take(assembly_sample_count)
                        .collect_vec();
                };
            }
        }
    }

    fn bin_contigs(&self, output: &str, m: &clap::ArgMatches) {
        match self {
            VariantMatrix::VariantContigMatrix { .. } => {
                external_command_checker::check_for_flock();
                let n_components = m.value_of("n-components").unwrap().parse().unwrap();
                // let min_cluster_size

                let cmd_string = format!(
                    "flock bin \
                    --assembly {} \
                    --input {} \
                    --kmer_frequencies {} \
                    --variant_rates {} \
                    --min_contig_size {} \
                    --n_neighbors {} \
                    --n_components {} \
                    --min_samples {} \
                    --min_dist 0 \
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
                    match m.is_present("variant-rates") {
                        true => m.value_of("variant-rates").unwrap().to_string(),
                        false => format!("{}/rosella_variant_rates.tsv", &output),
                    },
                    m.value_of("min-contig-size").unwrap(),
                    m.value_of("n-neighbors").unwrap(),
                    std::cmp::max(n_components, 2),
                    m.value_of("min-samples").unwrap(),
                    format!("{}/", &output),
                    m.value_of("threads").unwrap(),
                );

                command::finish_command_safely(
                    std::process::Command::new("bash")
                        .arg("-c")
                        .arg(&cmd_string)
                        .stderr(std::process::Stdio::piped())
                        // .stdout(std::process::Stdio::piped())
                        .spawn()
                        .expect("Unable to execute bash"),
                    "flock",
                );
            }
        }
    }

    fn finalize_bins(&mut self, output: &str, m: &clap::ArgMatches) {
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                target_lengths,
                sample_names,
                coverages,
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
                        if target_lengths.get(tid).unwrap() < &min_contig_size {
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

                // remove bins
                for bin in removed_bins {
                    bins.remove(&bin).unwrap();
                }

                // If we have enough samples, we can attempt to rescue the contigs we just
                // removed.
                // If n_samples is < 3, then flock will handle the recollection of large contigs
                // using soft clutering
                if sample_names.len() >= 3 {
                    // Iterate through removed contigs and check their average correlation with
                    // other bins
                    correlate_with_bins(
                        &removed_contigs[..],
                        &mut bins,
                        &coverages,
                        sample_names.len(),
                    );

                    // Now we use a similar method to above to rescue small contigs
                    correlate_with_bins(
                        &small_contigs[..],
                        &mut bins,
                        &coverages,
                        sample_names.len(),
                    );
                }

                *final_bins = bins
            }
        }
    }

    fn write_bins(&self, output: &str, reference: &str) {
        match self {
            VariantMatrix::VariantContigMatrix {
                final_bins,
                target_names,
                ..
            } => {
                final_bins.iter().for_each(|(bin, contigs)| {
                    let mut reference_file = generate_faidx(reference);
                    let mut writer = bio::io::fasta::Writer::to_file(format!(
                        "{}/rosella_bin_{}.fna",
                        &output, bin
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
        kmer_path: Option<&str>,
        rates_path: Option<&str>,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut target_names,
                ref mut target_lengths,
                ref mut coverages,
                ref mut variances,
                ref mut variant_rates,
                ref mut kfrequencies,
                ref mut sample_names,
                ..
            } => {
                match coverage_path {
                    Some(coverage_str) => {
                        // Read the coverage information in
                        let coverage_file = std::fs::File::open(coverage_str).unwrap();

                        let coverage_buffer = std::io::BufReader::new(coverage_file);
                        for (tid, line) in coverage_buffer.lines().enumerate() {
                            let line = line.unwrap();
                            let values = line.split('\t').collect_vec();

                            if tid == 0 {
                                for value in values[3..].iter().step_by(2) {
                                    let value = value.replace(".bam", "");
                                    sample_names.push(value);
                                }
                            } else {
                                target_names
                                    .entry(tid as i32 - 1)
                                    .or_insert(String::from(values[0]));
                                target_lengths
                                    .entry(tid as i32 - 1)
                                    .or_insert(values[1].parse().unwrap());

                                for (idx, value) in values[3..].iter().enumerate() {
                                    if idx % 2 == 0 {
                                        let idx = idx / 2;
                                        let contig_coverages = coverages
                                            .entry(tid as i32 - 1)
                                            .or_insert(vec![0.; sample_names.len()]);
                                        contig_coverages[idx as usize] = value.parse().unwrap();
                                    } else {
                                        let idx = idx / 2; // Rust defaults to floor division, so this should yield correct index
                                        let contig_variances = variances
                                            .entry(tid as i32 - 1)
                                            .or_insert(vec![0.; sample_names.len()]);
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
                        for (tid, line) in kmer_buffer.lines().enumerate() {
                            let line = line.unwrap();
                            let values = line.split('\t').collect_vec();

                            if tid == 0 {
                                for value in values[2..].iter() {
                                    kfrequencies
                                        .entry(value.as_bytes().to_vec())
                                        .or_insert(vec![0; target_names.len()]);
                                }
                            } else {
                                for (value, (_kmer, kmer_vec)) in
                                    values[2..].iter().zip(kfrequencies.iter_mut())
                                {
                                    kmer_vec[tid - 1] = value.parse().unwrap();
                                }
                            }
                        }
                    }
                    None => {
                        // Nothing to read, move on
                    }
                }

                match rates_path {
                    Some(rates_str) => {
                        // Read the coverage information in
                        let rates_file = std::fs::File::open(rates_str).unwrap();

                        let rates_buffer = std::io::BufReader::new(rates_file);
                        for (tid, line) in rates_buffer.lines().enumerate() {
                            let line = line.unwrap();
                            let values = line.split('\t').collect_vec();

                            if tid == 0 {
                                // Don't do anything in this case
                            } else {
                                for (idx, value) in values[2..].iter().enumerate() {
                                    if idx % 2 == 0 {
                                        let idx = idx / 2;
                                        let contig_rates = variant_rates
                                            .entry(tid as i32 - 1)
                                            .or_insert(HashMap::new());
                                        let snv_rates = contig_rates
                                            .entry(Var::SNV)
                                            .or_insert(vec![0.; sample_names.len()]);
                                        snv_rates[idx as usize] = value.parse().unwrap();
                                    } else {
                                        let idx = idx / 2; // Rust defaults to floor division, so this should yield correct index
                                        let contig_rates = variant_rates
                                            .entry(tid as i32 - 1)
                                            .or_insert(HashMap::new());
                                        let sv_rates = contig_rates
                                            .entry(Var::SV)
                                            .or_insert(vec![0.; sample_names.len()]);
                                        sv_rates[idx as usize] = value.parse().unwrap();
                                    }
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

pub fn correlate_with_bins(
    current_contigs: &[i32],
    bins: &mut HashMap<usize, Vec<i32>>,
    coverages: &HashMap<i32, Vec<f64>>,
    n: usize,
) {
    let to_add: Vec<(usize, i32)> = current_contigs
        .par_iter()
        .filter_map(|tid| {
            let (avg_correlation_s, avg_correlation_r) = channel();
            bins.par_iter()
                .for_each_with(avg_correlation_s, |avg_s, (bin, contigs)| {
                    let mut workspace = vec![0.; 2 * n];
                    let mut corr_sum = 0.;
                    contigs.iter().for_each(|other_tid| {
                        // Get TID of current index by indexing into large contigs, then
                        // uset that tid to get the coverage values
                        let current_cov = coverages.get(tid).unwrap();
                        let other_cov = coverages.get(other_tid).unwrap();
                        let spear =
                            spearman(&current_cov[..], 1, &other_cov[..], 1, n, &mut workspace);
                        corr_sum += spear;
                    });
                    let avg_corr = corr_sum / contigs.len() as f64;
                    avg_s.send((*bin, avg_corr)).unwrap();
                });
            let avg_correlations: Vec<(usize, f64)> = avg_correlation_r.iter().collect_vec();
            let mut max_bin = 0;
            let mut max_corr = 0.;

            // Find the maximum bin correlation
            for (bin, corr) in avg_correlations.iter() {
                if corr > &max_corr {
                    max_corr = *corr;
                    max_bin = *bin;
                }
            }

            // If the correlation is sufficient, place that contig in with that bin
            if max_corr >= 0.8 {
                Some((max_bin, *tid))
            } else {
                None
            }
        })
        .collect();

    for (max_bin, tid) in to_add.iter() {
        let mut bin = bins.entry(*max_bin).or_insert(Vec::new());
        bin.push(*tid);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::variants;
    use std::collections::HashSet;

    fn create_base(
        ref_sequence: &Vec<u8>,
        var_char: u8,
        pos: i64,
        sample_count: usize,
        depth: Vec<i32>,
        reference: Vec<i32>,
    ) -> Base {
        let mut total_depth = Vec::with_capacity(depth.len());
        for (v, r) in depth.iter().zip(reference.iter()) {
            total_depth.push(*v + *r);
        }
        Base {
            tid: 0,
            pos,
            refr: ref_sequence[pos as usize..(pos as usize + 1)].to_vec(),
            variant: Variant::SNV(var_char),
            depth: depth.clone(),
            truedepth: depth,
            totaldepth: total_depth,
            quals: vec![0.; sample_count],
            referencedepth: reference,
        }
    }

    #[test]
    fn test_stats_and_mat() {
        let ref_sequence = "ATGAAACCCGGGTTTTAA".as_bytes().to_vec();
        let sample_count = 2;

        let mut variant_abundances: HashMap<i64, HashMap<Variant, Base>> = HashMap::new();

        // Setup variant matrix
        let mut var_mat = VariantMatrix::new_matrix(2);

        // Create fake variants
        let var_1 = create_base(
            &ref_sequence,
            "G".bytes().nth(0).unwrap(),
            7,
            2,
            vec![0, 5],
            vec![5, 0],
        );
        let var_2 = create_base(
            &ref_sequence,
            "C".bytes().nth(0).unwrap(),
            11,
            2,
            vec![0, 5],
            vec![5, 0],
        );
        let var_3 = create_base(
            &ref_sequence,
            "A".bytes().nth(0).unwrap(),
            13,
            2,
            vec![0, 5],
            vec![5, 0],
        );
        let var_4 = create_base(
            &ref_sequence,
            "C".bytes().nth(0).unwrap(),
            14,
            2,
            vec![0, 5],
            vec![5, 0],
        );

        let mut ups_and_downs = vec![0; ref_sequence.len()];
        ups_and_downs[0] = 5;
        ups_and_downs[ref_sequence.len() - 1] = -5;

        // Setup sample 1
        let mut var_stats = VariantStats::new_contig_stats(0., 5., 0);

        var_stats.add_contig(
            Some(&mut variant_abundances),
            0,
            0,
            b"test".to_vec(),
            ref_sequence.len(),
            0,
            vec![10., 10., 0.],
            ups_and_downs,
        );
        var_mat.add_contig(var_stats, 2, 0);

        {
            // Add variants in
            let hash = variant_abundances.entry(7).or_insert(HashMap::new());
            hash.insert(var_1.variant.clone(), var_1);

            let hash = variant_abundances.entry(11).or_insert(HashMap::new());
            hash.insert(var_2.variant.clone(), var_2);

            let hash = variant_abundances.entry(13).or_insert(HashMap::new());
            hash.insert(var_3.variant.clone(), var_3);

            let hash = variant_abundances.entry(14).or_insert(HashMap::new());
            hash.insert(var_4.variant.clone(), var_4);
        }
        // Add sample 2
        let mut ups_and_downs = vec![0; ref_sequence.len()];
        ups_and_downs[0] = 5;
        ups_and_downs[ref_sequence.len() - 1] = -5;
        let mut var_stats = VariantStats::new_contig_stats(0., 5., 0);
        var_stats.add_contig(
            Some(&mut variant_abundances),
            0,
            0,
            b"test".to_vec(),
            ref_sequence.len(),
            1,
            vec![10., 10., 0.],
            ups_and_downs,
        );

        var_mat.add_contig(var_stats, 2, 1);

        let mut ref_map = HashMap::new();
        ref_map.insert(0, "test".to_string());
    }

    #[test]
    fn test_variant_removal() {
        let ref_sequence = "ATGAAACCCGGGTTTTAA".as_bytes().to_vec();

        // Setup variant matrix
        let mut var_mat = VariantMatrix::new_matrix(1);
        var_mat.add_info(0, "o".as_bytes().to_vec(), 20);
        // Create fake variants
        let var_1 = create_base(
            &ref_sequence,
            "G".bytes().nth(0).unwrap(),
            7,
            2,
            vec![5],
            vec![5],
        );
        let var_2 = create_base(
            &ref_sequence,
            "C".bytes().nth(0).unwrap(),
            11,
            2,
            vec![5],
            vec![5],
        );
        let var_3 = create_base(
            &ref_sequence,
            "A".bytes().nth(0).unwrap(),
            13,
            2,
            vec![5],
            vec![5],
        );
        let var_4 = create_base(
            &ref_sequence,
            "C".bytes().nth(0).unwrap(),
            14,
            2,
            vec![30],
            vec![30],
        );

        {
            // Add variants in
            var_mat.add_variant_to_matrix(0, &var_1, 0);

            var_mat.add_variant_to_matrix(0, &var_2, 0);

            var_mat.add_variant_to_matrix(0, &var_3, 0);

            var_mat.add_variant_to_matrix(0, &var_4, 0);
        }
        var_mat.remove_variants(0, 0, vec![30., 30.]);

        // retrieve variants
        let variants = var_mat.variants_of_contig(0).unwrap();
        assert_eq!(variants.len(), 4);

        var_mat.remove_variants(0, 0, vec![5., 5.]);

        // retrieve variants
        let variants = var_mat.variants_of_contig(0).unwrap();
        assert_eq!(variants.len(), 3);
    }
}
