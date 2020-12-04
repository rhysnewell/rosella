use rayon::prelude::*;
use std::collections::HashMap;
use std::io::prelude::*;
use std::str;

use std::fs::OpenOptions;
use std::path::Path;

use model::variants::*;

pub enum VariantStats {
    VariantContigStats {
        variants: HashMap<i64, HashMap<Variant, Base>>,
        variant_count: Vec<f64>,
        depth: Vec<i32>,
        tid: i32,
        total_indels: usize,
        target_name: Vec<u8>,
        target_len: f64,
        coverage: f64,
        variance: f64,
        observed_contig_length: u32,
        num_covered_bases: i32,
        num_mapped_reads: u64,
        total_mismatches: u64,
        contig_end_exclusion: u64,
        min: f64,
        max: f64,
        method: String,
        regression: (f64, f64, f64),
    },
}

#[allow(unused)]
impl VariantStats {
    pub fn new_contig_stats(min: f64, max: f64, contig_end_exclusion: u64) -> VariantStats {
        VariantStats::VariantContigStats {
            variants: HashMap::new(),
            tid: 0,
            variant_count: vec![],
            depth: vec![],
            total_indels: 0,
            target_name: vec![],
            target_len: 0.0,
            coverage: 0.00,
            variance: 0.00,
            observed_contig_length: 0,
            num_covered_bases: 0,
            contig_end_exclusion: contig_end_exclusion,
            num_mapped_reads: 0,
            total_mismatches: 0,
            min: min,
            max: max,
            method: "".to_string(),
            regression: (0., 0., 0.),
        }
    }
}

pub trait VariantFunctions {
    fn setup(&mut self);

    fn len(&mut self) -> usize;

    fn add_contig(
        &mut self,
        variant_map: Option<&mut HashMap<i64, HashMap<Variant, Base>>>,
        tid: i32,
        total_indels_in_contig: usize,
        contig_name: Vec<u8>,
        contig_len: usize,
        sample_idx: usize,
        coverages: Vec<f64>,
        ups_and_downs: Vec<i32>,
    );

    //    /// Prints out variant info for current contig
    //    fn print_variants(&mut self, ref_sequence: &Vec<u8>, stoit_name: &str);
}

impl VariantFunctions for VariantStats {
    fn setup(&mut self) {
        match self {
            VariantStats::VariantContigStats {
                ref mut variants,
                ref mut variant_count,
                ref mut depth,
                ref mut tid,
                ref mut total_indels,
                ref mut target_name,
                ref mut target_len,
                ref mut coverage,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            } => {
                *variants = HashMap::new();
                *variant_count = Vec::new();
                *depth = vec![];
                *tid = 0;
                *total_indels = 0;
                *target_name = vec![];
                *target_len = 0.0;
                *coverage = 0.00;
                *num_covered_bases = 0;
                *num_mapped_reads = 0;
            }
        }
    }

    #[allow(unused)]
    fn len(&mut self) -> usize {
        match self {
            VariantStats::VariantContigStats {
                ref mut variants, ..
            } => variants.len(),
        }
    }

    #[allow(unused)]
    fn add_contig(
        &mut self,
        variant_map: Option<&mut HashMap<i64, HashMap<Variant, Base>>>,
        target_id: i32,
        total_indels_in_contig: usize,
        contig_name: Vec<u8>,
        contig_len: usize,
        sample_idx: usize,
        coverages: Vec<f64>,
        ups_and_downs: Vec<i32>,
    ) {
        match self {
            VariantStats::VariantContigStats {
                ref mut variants,
                ref mut variant_count,
                ref mut depth,
                ref mut tid,
                ref mut total_indels,
                ref mut target_name,
                ref mut target_len,
                ref mut variance,
                ref mut coverage,
                ref mut method,
                ..
            } => {
                *tid = target_id;
                *total_indels = total_indels_in_contig;
                *target_name = contig_name;
                *target_len = contig_len as f64;
                *coverage = coverages[1] as f64;
                *variance = coverages[2] as f64;
                *method = method.to_string();
                *variants = match variant_map {
                    Some(map) => map.clone(),
                    _ => HashMap::new(),
                };

                // Cumulative sum of ups and downs vec to get depth
                *depth = ups_and_downs
                    .iter()
                    .scan(0, |acc, &x| {
                        *acc = *acc + x;
                        Some(*acc)
                    })
                    .collect();

                debug!(
                    "new contig added {} with coverage {} and variance {}",
                    tid, coverage, variance
                );
            }
        }
    }
}

// helper function to get the index of condensed matrix from it square form
#[allow(dead_code)]
fn condensed_index(i: usize, j: usize, n: usize) -> Option<usize> {
    if i == j {
        return None;
    } else {
        return Some(n * i - i * (i + 1) / 2 + j - 1 - i);
        //        return Some(n*(n-1)/2 - (n - row_i)*(n - row_i - 1)/2 + col_j - row_i - 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use coverm::bam_generator::*;
    use coverm::genome_exclusion::*;
    use coverm::mapping_parameters::*;
    use coverm::shard_bam_reader::*;
    use std::fs::File;
    use std::str;

    //    fn test_with_stream<R: NamedBamReader + Send,
    //        G: NamedBamReaderGenerator<R> + Send>(
    //        expected: &str,
    //        bam_readers: Vec<G>,
    //        mut reference: bio::io::fasta::IndexedReader<File>,
    //        proper_pairs_only: bool,
    //        n_threads: usize,
    //        coverage_fold: f32,
    //        min_var_depth: usize,
    //        min: f32,
    //        max: f32,
    //        mode: &str,
    //        include_indels: bool,
    //        include_soft_clipping: bool) {
    ////        let mut stream = Cursor::new(Vec::new());
    //        {
    //            reads_mapped_vec = pileup_variants(
    //                bam_readers,
    //                &mut coverage_taker,
    //                coverage_estimators,
    //                print_zero_coverage_contigs,
    //                flag_filters,
    //                false,
    //                );
    //        }
    ////        assert_eq!(expected, str::from_utf8(stream.get_ref()).unwrap());
    //    }
}
