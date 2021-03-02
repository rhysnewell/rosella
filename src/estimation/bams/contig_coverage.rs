use bio::io::fasta::IndexedReader;
use coverm::bam_generator::*;
use estimation::variant_matrix::*;
use rayon::prelude::*;
use rust_htslib::bam::{self, record::Cigar};
use std::fs::File;
use utils::*;

use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::FlagFilter;
use std::str;

/// Process all reads in a BAM file
#[allow(unused)]
pub fn contig_coverage<'b, R: IndexedNamedBamReader>(
    mut bam_generated: R,
    sample_idx: usize,
    sample_count: usize,
    coverage_estimators: &'b mut Vec<CoverageEstimator>,
    variant_matrix: &'b mut VariantMatrix,
    split_threads: usize,
    m: &'b clap::ArgMatches,
    output_prefix: &str,
    coverage_fold: f32,
    min_var_depth: usize,
    contig_end_exclusion: u64,
    min: f32,
    max: f32,
    mode: &'b str,
    include_soft_clipping: bool,
    include_indels: bool,
    flag_filters: &'b FlagFilter,
    mapq_threshold: u8,
    method: &'b str,
    readtype: ReadType,
    reference_file: &mut IndexedReader<File>,
    tid: u32,
) {
    // Sample name from bam file name
    let stoit_name = bam_generated.name().to_string().replace("/", ".");
    variant_matrix.add_sample_name(stoit_name.to_string(), sample_idx, readtype);

    // Set bam reading threads
    bam_generated.set_threads(split_threads);
    debug!("Managed to set threads.");

    // Retrieve bam header struct: Contains info about contigs and mappings
    let header = bam_generated.header().clone(); // bam header
    let target_name = header.target_names()[tid as usize]; // contig names

    let mut record: bam::record::Record = bam::Record::new(); // Empty bam record
    let mut ups_and_downs: Vec<i32> = Vec::new(); // Mosdepth coverage array
    let mut contig_name = Vec::new(); // Current contig name
    let mut contig_name_str = String::new(); // Human readable
    let mut num_mapped_reads_total: u64 = 0; // Reads mapped for this sample
    let mut num_mapped_reads_in_current_contig: u64 = 0; // Reads mapped in this sample for current contig
    let mut total_edit_distance_in_current_contig: u64 = 0; // Edit distance calculated from CIGAR strings

    let mut ref_seq: Vec<u8> = Vec::new(); // container for reference contig
    let mut last_tid: i32 = -2; // no such tid in a real BAM file
    let mut total_indels_in_current_contig = 0;

    let mut skipped_reads = 0; // for record in records

    // for each genomic position, retrive info for CoverM depth calculations.
    // Check for variant at that location. Check if mapping supports that variant.
    // Retrieve read ids for each variant
    // for (tid, target) in target_names.iter().enumerate() {
    let target_name = String::from_utf8(target_name.to_vec()).unwrap();
    // Contig names have reference appended to start
    let target_len = header.target_len(tid).unwrap(); // contig length
    variant_matrix.add_info(tid as usize, target_name.as_bytes().to_vec(), target_len);

    let mut ups_and_downs: Vec<i32> = vec![0; target_len as usize]; // Populate mosdepth array to current contig length

    bam_generated.fetch((tid)); // Retrieve current contig from BAM

    while bam_generated.read(&mut record) == true
    // Read records into the empty recorde
    {
        if (!flag_filters.include_supplementary
            && record.is_supplementary()
            && readtype != ReadType::Long) // We want supp alignments for longreads
            || (!flag_filters.include_secondary
            && record.is_secondary())
            || (!flag_filters.include_improper_pairs
            && !record.is_proper_pair()
            && readtype != ReadType::Long)
        // Check against filter flags and current sample type
        {
            skipped_reads += 1;
            continue;
        }

        // if reference has changed, print the last record
        let tid = record.tid();
        if !record.is_unmapped() {
            // if mapped
            if record.seq().len() == 0 {
                skipped_reads += 1;
                continue;
            } else if record.mapq() < mapq_threshold {
                skipped_reads += 1;
                continue;
            }

            if !record.is_supplementary() {
                num_mapped_reads_in_current_contig += 1;
            }

            // for each chunk of the cigar string
            let mut cursor: usize = record.pos() as usize;
            let read_len = record.seq().len();
            for cig in record.cigar().iter() {
                match cig {
                    Cigar::Match(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                        // if M, X, or = increment start and decrement end index
                        ups_and_downs[cursor] += 1;
                        let final_pos = cursor + cig.len() as usize;

                        if final_pos < ups_and_downs.len() {
                            // True unless the read hits the contig end.
                            ups_and_downs[final_pos] -= 1;
                        }
                        cursor += cig.len() as usize;
                    }
                    Cigar::Del(del) => {
                        cursor += cig.len() as usize;
                        total_indels_in_current_contig += cig.len() as u64;
                    }
                    Cigar::RefSkip(del) => {
                        cursor += cig.len() as usize;
                    }
                    Cigar::Ins(ins) => {
                        total_indels_in_current_contig += cig.len() as u64;
                    }
                    Cigar::SoftClip(_) => {}
                    Cigar::HardClip(_) | Cigar::Pad(_) => {}
                }
            }
            // Determine the number of mismatching bases in this read by
            // looking at the NM tag.
            total_edit_distance_in_current_contig += match record.aux("NM".as_bytes()) {
                Some(aux) => aux.integer() as u64,
                None => {
                    panic!(
                        "Mapping record encountered that does not have an 'NM' \
                    auxiliary tag in the SAM/BAM format. This is required \
                    to work out some coverage statistics"
                    );
                }
            };
        }
    }
    let contig_len = header.target_len(tid as u32).expect("Corrupt BAM file?") as usize;
    contig_name = target_name.as_bytes().to_vec();

    let total_mismatches = total_edit_distance_in_current_contig - total_indels_in_current_contig;

    process_previous_contigs_var(
        mode,
        readtype,
        tid as i32,
        ups_and_downs,
        coverage_estimators,
        variant_matrix,
        sample_idx,
        total_mismatches,
        num_mapped_reads_in_current_contig,
        sample_count,
    );

    num_mapped_reads_total += num_mapped_reads_in_current_contig;

    // remove tmp file name from sample id
    let stoit_name = match stoit_name.to_string().contains(".tmp") {
        true => &stoit_name[15..],
        false => &stoit_name,
    };

    debug!(
        "In sample '{}', found {} reads mapped out of {} total ({:.*}%)", // and filtered {}",
        stoit_name,
        num_mapped_reads_total,
        bam_generated.num_detected_primary_alignments(),
        2,
        (num_mapped_reads_total * 100) as f64
            / bam_generated.num_detected_primary_alignments() as f64,
        // skipped_reads
    );

    bam_generated.finish();
}

#[allow(unused)]
pub fn process_previous_contigs_var(
    mode: &str,
    read_type: ReadType,
    last_tid: i32,
    ups_and_downs: Vec<i32>,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    variant_matrix: &mut VariantMatrix,
    sample_idx: usize,
    total_mismatches: u64,
    num_mapped_reads_in_current_contig: u64,
    sample_count: usize,
) {
    if last_tid != -2 {
        coverage_estimators
            .par_iter_mut()
            .for_each(|estimator| estimator.setup());

        coverage_estimators.par_iter_mut().for_each(|estimator| {
            estimator.add_contig(
                &ups_and_downs,
                num_mapped_reads_in_current_contig,
                total_mismatches,
            )
        });

        let coverages: Vec<f64> = coverage_estimators
            .iter_mut()
            .map(|estimator| estimator.calculate_coverage(&vec![0]) as f64)
            .collect();

        match mode {
            "bin" => {
                // Add samples contig information to main struct
                debug!("Adding in new info for contig...");
                variant_matrix.add_contig(last_tid, coverages, sample_count, sample_idx, read_type);
            }
            _ => {
                panic!("unknown mode {}", mode);
            }
        }
    }
}
