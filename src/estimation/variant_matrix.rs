use bio::io::gff;
use bio::stats::{
    bayesian,
    probs::{LogProb, PHREDProb, Prob},
};
use bird_tool_utils::command;
use estimation::contig_variants::*;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use itertools::izip;
use itertools::Itertools;
use model::variants::*;
use ndarray::prelude::*;
use ndarray_npy::{read_npy, write_npy};
use ordered_float::NotNan;
use rayon::current_num_threads;
use rayon::prelude::*;
use rust_htslib::bcf::{self};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::process::Stdio;
use std::str;
use std::sync::{Arc, Mutex};
use tempfile;
use utils::generate_faidx;

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
    fn remove_variants(
        &mut self,
        sample_idx: usize,
        contig_info: HashMap<usize, Vec<f64>>,
    );

    /// Per sample FDR calculation
    fn remove_false_discoveries(&mut self, alpha: f64, reference: &str);

    /// Returns the alleles at the current position
    /// as a mutable reference
    fn variants(
        &mut self,
        tid: i32,
        pos: i64,
    ) -> Option<&mut HashMap<Variant, Base>>;

    /// Returns all variants found within a contig
    fn variants_of_contig(
        &mut self,
        tid: i32,
    ) -> Option<&mut HashMap<i64, HashMap<Variant, Base>>>;

    /// Takes [VariantStats](contig_variants/VariantStats) struct for single contig and adds to
    /// [VariantMatrix](VariantMatrix)
    fn add_contig(
        &mut self,
        variant_stats: VariantStats,
        sample_count: usize,
        sample_idx: usize,
    );

    /// Converts all variants into variant::Var format
    fn generate_distances(&mut self);

    /// Prints the per reference and per contig variant information e.g. How many SNPs were seen
    /// along a contig over the given window size on average
    fn print_variant_stats(
        &self,
        window_size: f64,
        output_prefix: &str,
    );

    fn calculate_sample_distances(
        &self,
        output_prefix: &str,
    );

    fn write_vcf(&self, output_prefix: &str,);
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
            VariantMatrix::VariantContigMatrix {
                variant_counts,
                ..
            } => {
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
                match &sample_names[sample_idx][..4] {
                    ".tmp" => &sample_names[sample_idx][15..],
                    _ => &sample_names[sample_idx],
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

    fn add_variant_to_matrix(
        &mut self,
        sample_idx: usize,
        base: &Base,
        tid: usize,
    ) {
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
                let contig_variants = all_variants
                    .entry(tid as i32)
                    .or_insert(HashMap::new());

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

    fn remove_variants(
        &mut self,
        sample_idx: usize,
        contig_info: HashMap<usize, Vec<f64>>,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut all_variants,
                ..
            } => {

                for (tid, stats) in contig_info.iter() {
                    let upper_limit = stats[0] + stats[1];
                    let lower_limit = stats[0] - stats[1];
                    match all_variants.get_mut(&(*tid as i32)) {
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
                                                    let window_variants = contig_variants
                                                        .get(&window)
                                                        .unwrap();
                                                    for (window_variant, _) in
                                                        window_variants.iter()
                                                    {
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
                                    if (total_depth < lower_limit
                                        || total_depth > upper_limit)
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
    }

    fn remove_false_discoveries(&mut self, alpha: f64, reference: &str) {
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
                                    NotNan::new(allele_probs)
                                        .expect("Unable to convert to NotNan"),
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

    fn variants(
        &mut self,
        tid: i32,
        pos: i64,
    ) -> Option<&mut HashMap<Variant, Base>> {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut all_variants,
                ..
            } => {
                match all_variants.get_mut(&tid) {
                    Some(contig_variants) => contig_variants.get_mut(&pos),
                    None => None,
                }
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
            } => {
                all_variants.get_mut(&tid)
            },
        }
    }

    fn add_contig(
        &mut self,
        variant_stats: VariantStats,
        sample_count: usize,
        sample_idx: usize,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                ref mut coverages,
                ref mut all_variants,
                target_names,
                //                ref mut target_lengths,
                ref mut variances,
                ref mut depths,
                ..
            } => {
                match variant_stats {
                    VariantStats::VariantContigStats {
                        tid,
                        //                        target_name,
                        //                        target_len,
                        coverage,
                        variance,
                        depth,
                        //                        variations_per_n,
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
                }
            }
        }
    }

    fn generate_distances(&mut self) {
        match self {
            VariantMatrix::VariantContigMatrix {
                all_variants,
                sample_names,
                coverages,
                ref mut variant_info,
                ref mut geom_mean_var,
                ref mut geom_mean_dep,
                ref mut geom_mean_frq,
                ..
            } => {
                let sample_count = sample_names.len();


                // A vector the length of the number of samples
                // Cumulatively calculates the product of variant depths

                let variant_info_vec = Arc::new(Mutex::new(Vec::new()));

                // A vector the length of the number of samples
                // Cumulatively calculates the product of variant depths
                let geom_mean_v = Arc::new(Mutex::new(vec![1. as f64; sample_count as usize]));

                // product of total depth
                let geom_mean_d = Arc::new(Mutex::new(vec![1. as f64; sample_count as usize]));

                // product of reference frequency
                let geom_mean_f = Arc::new(Mutex::new(vec![1. as f64; sample_count as usize]));

                // get basic variant info and store as variant::Var

                all_variants
                    .iter_mut()
                    .for_each(|(tid, variant_abundances)| {
                        variant_abundances
                            .par_iter_mut()
                            .for_each(|(position, hash)| {
                                // loop through each position that has variants ignoring positions that
                                // only contained the reference in all samples
                                if hash.keys().len() > 0 {
                                    for (variant, base_info) in hash.iter_mut() {
                                        match variant {
                                            Variant::SNV(_) | Variant::None => {
                                                let _abundance: f64 = 0.;
                                                let _mean_var: f64 = 0.;
                                                let mut rel_abund =
                                                    vec![0.0; sample_count as usize];

                                                let depth_sum: i32 =
                                                    base_info.truedepth.iter().sum();
                                                if depth_sum > 0 {
                                                    // Get the mean abundance across samples
                                                    for index in (0..sample_count).into_iter() {
                                                        let mut geom_mean_v =
                                                            geom_mean_v.lock().unwrap();
                                                        let mut geom_mean_d =
                                                            geom_mean_d.lock().unwrap();
                                                        let mut geom_mean_f =
                                                            geom_mean_f.lock().unwrap();

                                                        let mut var_depth =
                                                            base_info.truedepth[index] as f64;

                                                        let total_depth =
                                                            base_info.totaldepth[index] as f64;
                                                        //                                                println!("var_depth {} tot {}", var_depth, total_depth);
                                                        //                                                base_info.freq[index] = ;
                                                        if total_depth <= 0. {
                                                            rel_abund[index] = var_depth / 1.;
                                                        } else {
                                                            rel_abund[index] =
                                                                var_depth / total_depth;
                                                        }

                                                        geom_mean_v[index] +=
                                                            (var_depth + 1.).ln();
                                                        geom_mean_d[index] +=
                                                            (total_depth + 1.).ln();
                                                        geom_mean_f[index] += ((var_depth
                                                            + 1.)
                                                            / (total_depth + 1.))
                                                            .ln();
                                                    }

                                                    let mut variant_info_vec =
                                                        variant_info_vec.lock().unwrap();
                                                    //                                            base_info.rel_abunds = rel_abund;
                                                    let point = Var {
                                                        pos: *position,
                                                        var: variant.clone(),
                                                        deps: base_info.totaldepth.clone(),
                                                        vars: base_info.truedepth.clone(),
                                                        //                                                rel_abunds: rel_abund,
                                                        tid: *tid,
                                                        reads: base_info.reads.clone(),
                                                    };
                                                    variant_info_vec.push(point);
                                                }
                                            }
                                            _ => {
                                                let mut rel_abund =
                                                    vec![0.0; sample_count as usize];
                                                let depth_sum: i32 =
                                                    base_info.truedepth.iter().sum();
                                                if depth_sum > 0 {
                                                    // Get the mean abundance across samples
                                                    for index in (0..sample_count).into_iter() {
                                                        let mut geom_mean_v =
                                                            geom_mean_v.lock().unwrap();
                                                        let mut geom_mean_d =
                                                            geom_mean_d.lock().unwrap();
                                                        let mut geom_mean_f =
                                                            geom_mean_f.lock().unwrap();

                                                        let mut var_depth =
                                                            base_info.depth[index] as f64;
                                                        if var_depth <= 0. {
                                                            var_depth = base_info.truedepth
                                                                [index]
                                                                as f64;
                                                            base_info.depth[index] =
                                                                base_info.truedepth[index];
                                                        }
                                                        let total_depth =
                                                            base_info.totaldepth[index] as f64;
                                                        //                                                base_info.freq[index] = ;
                                                        if total_depth <= 0. {
                                                            rel_abund[index] = var_depth / (1.);
                                                        } else {
                                                            rel_abund[index] =
                                                                var_depth / total_depth;
                                                        }

                                                        geom_mean_v[index] +=
                                                            (var_depth + 1.).ln();
                                                        geom_mean_d[index] +=
                                                            (total_depth + 1.).ln();
                                                        geom_mean_f[index] += ((var_depth
                                                            + 1.)
                                                            / (total_depth + 1.))
                                                            .ln();
                                                    }

                                                    let mut variant_info_vec =
                                                        variant_info_vec.lock().unwrap();
                                                    //                                            base_info.rel_abunds = rel_abund;
                                                    let point = Var {
                                                        pos: *position,
                                                        var: variant.clone(),
                                                        deps: base_info.totaldepth.clone(),
                                                        vars: base_info.truedepth.clone(),
                                                        //                                                rel_abunds: rel_abund,
                                                        tid: *tid,
                                                        reads: base_info.reads.clone(),
                                                    };
                                                    variant_info_vec.push(point);
                                                }
                                            }
                                        }
                                    }
                                }
                            });
                    });

                    // Add variant info for reference to hashmap
                    let variant_info_vec = variant_info_vec.lock().unwrap().clone();

                    // Helper fn to calculate geom mean from sum of log values
                    let geom_mean = |input: &Vec<f64>| -> Vec<f64> {
                        let output = input
                            .iter()
                            .map(|sum| (sum / variant_info_vec.len() as f64).exp())
                            .collect::<Vec<f64>>();
                        return output;
                    };

                    debug!(
                        "geoms {:?} {:?} {:?}",
                        geom_mean_d, geom_mean_v, geom_mean_f
                    );

                    let geom_mean_v = geom_mean_v.lock().unwrap().clone();
                    let geom_mean_v = geom_mean(&geom_mean_v);
                    debug!("Geom Mean Var {:?}", geom_mean_v);

                    let geom_mean_d = geom_mean_d.lock().unwrap().clone();
                    let geom_mean_d = geom_mean(&geom_mean_d);
                    debug!("Geom Mean Dep {:?}", geom_mean_d);

                    let geom_mean_f = geom_mean_f.lock().unwrap().clone();
                    let geom_mean_f = geom_mean(&geom_mean_f);
                    debug!("Geom Mean Frq {:?}", geom_mean_f);
                    debug!(
                        "geoms {:?} {:?} {:?}",
                        geom_mean_d, geom_mean_v, geom_mean_f
                    );


                *variant_info = variant_info_vec;
                *geom_mean_var = geom_mean_v;
                *geom_mean_dep = geom_mean_d;
                *geom_mean_frq = geom_mean_f;
            }
        }
    }

    fn print_variant_stats(
        &self,
        window_size: f64,
        output_prefix: &str,
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                target_lengths,
                sample_names,
                all_variants,
                variant_sums,
                variant_info,
                ..
            } => {
                let file_name = format!("{}/rosella_summary.tsv", &output_prefix);
                let snp_locs =
                    format!("{}/rosella_snp_locations.tsv", &output_prefix,);

                let file_path = Path::new(&file_name);

                let mut file_open = match File::create(file_path) {
                    Ok(fasta) => fasta,
                    Err(e) => {
                        println!("Cannot create file {:?}", e);
                        std::process::exit(1)
                    }
                };

                let snp_loc_path = Path::new(&snp_locs);
                let mut snp_loc_open = match File::create(snp_loc_path) {
                    Ok(tsv) => tsv,
                    Err(e) => {
                        println!("Cannot create file {:?}", e);
                        std::process::exit(1)
                    }
                };

                // Snp density summary start
                write!(file_open, "contigName\tcontigLen").unwrap();

                // Snp location start
                write!(snp_loc_open, "SNP\tchr\tpos").unwrap();
                debug!("Sample Names {:?}", sample_names);
                for sample_name in sample_names.iter() {
                    write!(
                        file_open,
                        "\t{}.snvsPer{}kb\t{}.svsPer{}kb\t{}.snvCount\t{}.svCount",
                        &sample_name,
                        &window_size,
                        &sample_name,
                        &window_size,
                        &sample_name,
                        &sample_name
                    )
                    .unwrap();
                }
                write!(file_open, "\n").unwrap();
                write!(snp_loc_open, "\n").unwrap();

                let mut snps = 0;
                for (tid, contig_name) in target_names.iter() {
                    let contig_len = target_lengths[&tid];
                    write!(file_open, "{}\t{}", contig_name, contig_len).unwrap();
                    match all_variants.get(tid) {
                        Some(variants_in_contig) => {
                            // Set up channels that receive a vector of values
                            // for each sample
                            let mut snps_cnt_vec = vec![0; sample_names.len()];
                            let svs_cnt_vec =
                                Arc::new(Mutex::new(vec![0; sample_names.len()]));

                            let window = contig_len / window_size;
                            variants_in_contig.iter().enumerate()
                                .for_each(
                                    |(index,
                                         (position, variants))| {
                                    // Get how many alleles are present at loci
                                    let alleles = variants.len();
                                    for (var, base) in variants {
                                        match var {
                                            Variant::SNV(_) => {
                                                write!(snp_loc_open, "SNP{}\t{}\t{}\n",
                                                       index, contig_name, position).unwrap();
                                                base.truedepth
                                                    .iter()
                                                    .enumerate()
                                                    .zip(base.depth.iter().enumerate())
                                                    .for_each(|((index, count_1), (_index_2, count_2))| {
                                                        if count_1 > &0 || count_2 > &0 {
                                                            snps_cnt_vec[index] += 1;
                                                        }
                                                    });
                                                snps += 1;
                                            },
                                            Variant::None => {
                                                // If biallelic or multiallelic then
                                                // I don't think we want the ref depth?
                                            },
                                            _ => {
                                                base.truedepth
                                                    .par_iter()
                                                    .enumerate()
                                                    .zip(base.depth.par_iter().enumerate())
                                                    .for_each(|((index, count_1), (_, count_2))| {
                                                        if count_1 > &0 || count_2 > &0 {
                                                            let mut svs_cnt_vec
                                                                = svs_cnt_vec.lock().unwrap();
                                                            svs_cnt_vec[index] += 1;
                                                        }
                                                    })
                                            }
                                        };
                                    };
                                });
                            let svs_cnt_vec = svs_cnt_vec.lock().unwrap().clone();
                            let snps_per_win: Vec<_> = snps_cnt_vec
                                .iter()
                                .map(|count| *count as f64 / window)
                                .collect();
                            let svs_per_win: Vec<_> = svs_cnt_vec
                                .iter()
                                .map(|count| *count as f64 / window)
                                .collect();

                            for (snp_w, svs_w, snp_c, svs_c) in izip!(
                                &snps_per_win,
                                &svs_per_win,
                                &snps_cnt_vec,
                                &svs_cnt_vec
                            ) {
                                write!(
                                    file_open,
                                    "\t{}\t{}\t{}\t{}",
                                    snp_w, svs_w, snp_c, svs_c
                                )
                                .unwrap();
                            }
                            write!(file_open, "\n").unwrap();
                        }
                            None => {
                                // Write out zeros for contigs with no variants
                                for (sample_idx, _sample_name) in
                                    sample_names.iter().enumerate()
                                {
                                    write!(
                                        file_open,
                                        "\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                        0., 0., 0., 0., 0., 0., 0.,
                                    )
                                    .unwrap();
                                }
                            }
                        }
                    }
            }
        }
    }

    fn calculate_sample_distances(
        &self,
        output_prefix: &str
    ) {
        match self {
            VariantMatrix::VariantContigMatrix {
                target_names,
                all_variants,
                sample_names,
                ..
            } => {
                let number_of_samples = sample_names.len();

                if all_variants.len() > 0 {
                    debug!(
                        "Generating variant adjacency matrix",
                    );

                    // The initialization of adjacency matrix n*n n=samples
                    let adjacency_matrix =
                        Arc::new(Mutex::new(vec![
                            vec![0; number_of_samples];
                            number_of_samples
                        ]));

                    let sample_totals = Arc::new(Mutex::new(vec![0; number_of_samples]));

                    debug!("Collecting variants...");

                    for (tid, contig_variants) in all_variants.iter() {
                        for (pos, pos_variants) in contig_variants.iter() {
                            if pos_variants.len() > 1 {
                                for (variant, base_info) in pos_variants.iter() {
                                    match variant {
                                        Variant::SNV(_)
                                        | Variant::MNV(_)
                                        | Variant::SV(_)
                                        | Variant::Insertion(_)
                                        | Variant::Inversion(_)
                                        | Variant::Deletion(_)
                                        | Variant::None => {
                                            // Update adjacency matrix for each variant
                                            // Evaluate each sample pairwise
                                            (0..number_of_samples)
                                                .into_iter()
                                                .combinations(2)
                                                .collect::<Vec<Vec<usize>>>()
                                                .into_par_iter()
                                                .for_each(|indices| {
                                                    let d1 =
                                                        base_info.truedepth[indices[0]];
                                                    let d2 =
                                                        base_info.truedepth[indices[1]];

                                                    if d1 > 0 && d2 > 0 {
                                                        // udpate both
                                                        let mut adjacency_matrix =
                                                            adjacency_matrix
                                                                .lock()
                                                                .unwrap();
                                                        let mut sample_totals =
                                                            sample_totals.lock().unwrap();
                                                        sample_totals[indices[0]] += 1;
                                                        sample_totals[indices[1]] += 1;

                                                        adjacency_matrix[indices[0]]
                                                            [indices[1]] += 1;
                                                        adjacency_matrix[indices[1]]
                                                            [indices[0]] += 1;
                                                    } else if d1 > 0 && d2 <= 0 {
                                                        // update only one
                                                        let mut sample_totals =
                                                            sample_totals.lock().unwrap();
                                                        sample_totals[indices[0]] += 1;
                                                    } else if d1 <= 0 && d2 > 0 {
                                                        let mut sample_totals =
                                                            sample_totals.lock().unwrap();
                                                        sample_totals[indices[1]] += 1;
                                                    }
                                                });
                                        }
                                        _ => {
                                            // do nothing
                                        }
                                    }
                                }
                            }
                        }
                    }

                    let adjacency_matrix = adjacency_matrix.lock().unwrap();
                    let sample_totals = sample_totals.lock().unwrap();

                    let file_name = format!(
                        "{}/rosella_adjacency_matrix.tsv",
                        &output_prefix,
                    );

                    let file_path = Path::new(&file_name);

                    let mut file_open = match File::create(file_path) {
                        Ok(fasta) => fasta,
                        Err(e) => {
                            println!("Cannot create file {:?}", e);
                            std::process::exit(1)
                        }
                    };

                    // Adjacency summary start
                    write!(file_open, "{: <20}", "").unwrap();
                    for sample_name in sample_names.iter() {
                        // remove tmp file name from sample id
                        let sample_name = match &sample_name[..4] {
                            ".tmp" => &sample_name[15..],
                            _ => &sample_name,
                        };
                        write!(file_open, "\t{: <20}", sample_name).unwrap();
                    }

                    write!(file_open, "\n").unwrap();

                    for (sample_idx_1, distance_vec) in adjacency_matrix.iter().enumerate()
                    {
                        let mut sample_name = &sample_names[sample_idx_1];
                        // remove tmp file name from sample id
                        let sample_name = match &sample_name[..4] {
                            ".tmp" => &sample_name[15..],
                            _ => &sample_name,
                        };
                        write!(file_open, "{: <20}", &sample_name,).unwrap();
                        for (sample_idx_2, count) in distance_vec.iter().enumerate() {
                            // Return the jaccard's similarity between sets of variants in
                            // samples

                            // Return shared variant count
                            let count = if sample_idx_1 == sample_idx_2 {
                                sample_totals[sample_idx_1] as f32
                                    / (number_of_samples - 1) as f32
                            } else {
                                *count as f32
                            };

                            write!(file_open, "\t{: <20}", count,).unwrap();
                        }
                        write!(file_open, "\n").unwrap();
                    }
                }
            }
        }
    }

    fn write_vcf(&self, output_prefix: &str) {
        match self {
            VariantMatrix::VariantContigMatrix {
                all_variants,
                target_names,
                target_lengths,
                sample_names,
                variant_info,
                ..
            } => {
                    if variant_info.len() > 0 {
                        // initiate header
                        let mut header = bcf::Header::new();
                        // Add program info
                        header.push_record(
                            format!("##source=lorikeet-v{}", env!("CARGO_PKG_VERSION")).as_bytes(),
                        );

                        debug!("samples {:?}", &sample_names);
                        for sample_name in sample_names.iter() {
                            // remove tmp file name from sample id
                            let sample_name = match &sample_name[..4] {
                                ".tmp" => &sample_name[15..],
                                _ => &sample_name,
                            };
                            header.push_sample(&sample_name.to_string().into_bytes()[..]);
                        }

                        // Add contig info
                        for (tid, contig_name) in target_names.iter() {
                            header.push_record(
                                format!(
                                    "##contig=<ID={}, length={}>",
                                    contig_name, target_lengths[&tid]
                                )
                                .as_bytes(),
                            );
                        }

                        // Add INFO flags
                        header.push_record(
                            format!(
                                "##INFO=<ID=TYPE,Number=A,Type=String,\
                    Description=\"The type of allele, either SNV, MNV, INS, DEL, or INV.\">"
                            )
                            .as_bytes(),
                        );

                        header.push_record(
                            format!(
                                "##INFO=<ID=SVLEN,Number=1,Type=Integer,\
                    Description=\"Length of structural variant\">"
                            )
                            .as_bytes(),
                        );

                        header.push_record(
                            b"##INFO=<ID=TDP,Number=A,Type=Integer,\
                            Description=\"Total observed sequencing depth\">",
                        );

                        header.push_record(
                            b"##INFO=<ID=TAD,Number=A,Type=Integer,\
                            Description=\"Total observed sequencing depth of alternative allele\">",
                        );

                        header.push_record(
                            b"##INFO=<ID=TRD,Number=A,Type=Integer,\
                            Description=\"Total observed sequencing depth of reference allele\">",
                        );

                        header.push_record(
                            b"##INFO=<ID=ST,Number=.,Type=Integer,\
                            Description=\"The strain IDs assigned to this variant\">",
                        );

                        // Add FORMAT flags
                        header.push_record(
                            format!(
                                "##FORMAT=<ID=DP,Number={},Type=Integer,\
                            Description=\"Observed sequencing depth in each sample\">",
                                sample_names.len()
                            )
                            .as_bytes(),
                        );

                        header.push_record(
                            format!("##FORMAT=<ID=AD,Number={},Type=Integer,\
                            Description=\"Observed sequencing depth of alternative allele in each sample\">",
                                    sample_names.len()).as_bytes(),
                        );

                        header.push_record(
                            format!("##FORMAT=<ID=RD,Number={},Type=Integer,\
                            Description=\"Observed sequencing depth of reference allele in each sample\">",
                                    sample_names.len()).as_bytes(),
                        );

                        header.push_record(
                            format!(
                                "##FORMAT=<ID=QA,Number={},Type=Float,\
                            Description=\"Quality scores for allele in each sample\">",
                                sample_names.len()
                            )
                            .as_bytes(),
                        );

                        let vcf_presort = tempfile::Builder::new()
                            .prefix("lorikeet-fifo")
                            .tempfile()
                            .expect("Unable to create VCF tempfile");

                        // Initiate writer
                        let mut bcf_writer = bcf::Writer::from_path(
                            format!(
                                "{}/rosella.vcf",
                                &output_prefix,
                            )
                            .as_str(),
                            &header,
                            true,
                            bcf::Format::VCF,
                        )
                        .expect(
                            format!("Unable to create VCF output: {}/rosella.vcf", output_prefix).as_str(),
                        );

                        bcf_writer.set_threads(current_num_threads()).unwrap();

                        for (tid, position_variants) in all_variants.iter() {
                            let contig_name = &target_names[&tid];
                            for (pos, variants) in position_variants.iter() {
                                for (variant, base) in variants.iter() {
                                    // Create empty record for this variant
                                    let mut record = bcf_writer.empty_record();
                                    record.set_rid(Some(
                                        bcf_writer
                                            .header()
                                            .name2rid(contig_name.as_bytes())
                                            .unwrap(),
                                    ));
                                    record.set_pos(*pos);

                                    // Sum the quality scores across samples
                                    let qual_sum = base.quals.par_iter().sum::<f64>();
                                    record.set_qual(qual_sum as f32);

                                    // Collect strain information
                                    let mut strains =
                                        base.genotypes.iter().cloned().collect::<Vec<i32>>();
                                    strains.sort();

                                    // Push info tags to record
                                    record
                                        .push_info_integer(b"TDP", &[base.totaldepth.iter().sum()]);
                                    record
                                        .push_info_integer(b"TAD", &[base.truedepth.iter().sum()]);
                                    record.push_info_integer(
                                        b"TRD",
                                        &[base.referencedepth.iter().sum()],
                                    );
                                    record.push_info_integer(b"ST", &strains[..]);

                                    // Push format flags to record
                                    record.push_format_integer(b"DP", &base.totaldepth[..]);
                                    record.push_format_integer(b"AD", &base.truedepth[..]);
                                    record.push_format_integer(b"RD", &base.referencedepth[..]);
                                    record.push_format_float(
                                        b"QA",
                                        &base.quals.iter().map(|p| *p as f32).collect_vec()[..],
                                    );

                                    let refr = &base.refr[..];

                                    match variant {
                                        Variant::SNV(alt) => {
                                            // Collect and set the alleles to record
                                            let alt = &[*alt];
                                            let mut collect_alleles = vec![refr, alt];
                                            record.set_alleles(&collect_alleles).unwrap();

                                            record.push_info_string(b"TYPE", &[b"SNV"]);

                                            bcf_writer
                                                .write(&record)
                                                .expect("Unable to write record");
                                        }
                                        Variant::MNV(alt) => {
                                            // Collect and set the alleles to record
                                            let alt = [&refr[..], &alt[..]].concat();
                                            let mut collect_alleles = vec![refr, &alt];
                                            record.set_alleles(&collect_alleles).unwrap();

                                            record.push_info_string(b"TYPE", &[b"MNV"]);
                                            record.push_info_integer(b"SVLEN", &[alt.len() as i32]);

                                            bcf_writer
                                                .write(&record)
                                                .expect("Unable to write record");
                                        }
                                        Variant::Inversion(alt) => {
                                            // Collect and set the alleles to record
                                            let alt = [&refr[..], &alt[..]].concat();
                                            let mut collect_alleles = vec![refr, &alt];
                                            record.set_alleles(&collect_alleles).unwrap();

                                            record.push_info_string(b"TYPE", &[b"INV"]);
                                            record.push_info_integer(b"SVLEN", &[alt.len() as i32]);

                                            bcf_writer
                                                .write(&record)
                                                .expect("Unable to write record");
                                        }
                                        Variant::Insertion(alt) => {
                                            // Collect and set the alleles to record
                                            let alt = [&refr[..], &alt[..]].concat();
                                            let mut collect_alleles = vec![refr, &alt];
                                            record.set_alleles(&collect_alleles).unwrap();

                                            record.push_info_string(b"TYPE", &[b"INS"]);
                                            record.push_info_integer(b"SVLEN", &[alt.len() as i32]);

                                            bcf_writer
                                                .write(&record)
                                                .expect("Unable to write record");
                                        }
                                        Variant::Deletion(alt) => {
                                            // Collect and set the alleles to record
                                            // Create deletion variant, attach refr to head of N
                                            record
                                                .set_alleles(&vec![
                                                    refr,
                                                    format!("{}N", refr[0] as char).as_bytes(),
                                                ])
                                                .unwrap();

                                            record.push_info_string(b"TYPE", &[b"DEL"]);
                                            record.push_info_integer(b"SVLEN", &[*alt as i32]);

                                            bcf_writer
                                                .write(&record)
                                                .expect("Unable to write record");
                                        }
                                        _ => {}
                                    }
                                }
                            }
                        }
                        debug!("Finished writing VCF file for {}", &output_prefix);
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
            genotypes: HashSet::new(),
            quals: vec![0.; sample_count],
            referencedepth: reference,
            freq: vec![0.; sample_count],
            rel_abunds: vec![0.; sample_count],
            reads: HashSet::new(),
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
        var_mat.add_contig(var_stats, 2, 0, 0);

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

        var_mat.generate_distances();
        let mut ref_map = HashMap::new();
        ref_map.insert(0, "test".to_string());
        let multi = Arc::new(MultiProgress::new());
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
            var_mat.add_variant_to_matrix(0, &var_1, 0, 0);

            var_mat.add_variant_to_matrix(0, &var_2, 0, 0);

            var_mat.add_variant_to_matrix(0, &var_3, 0, 0);

            var_mat.add_variant_to_matrix(0, &var_4, 0, 0);
        }
        let mut stats = HashMap::new();
        stats.insert(0, vec![30., 30.]);
        var_mat.remove_variants(0,  stats);

        // retrieve variants
        let variants = var_mat.variants_of_contig(0, 0).unwrap();
        assert_eq!(variants.len(), 4);

        let mut stats = HashMap::new();
        stats.insert(0, vec![5., 5.]);
        var_mat.remove_variants(0, stats);

        // retrieve variants
        let variants = var_mat.variants_of_contig(0, 0).unwrap();
        assert_eq!(variants.len(), 3);
    }
}
