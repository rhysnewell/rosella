use bio::io::fasta::IndexedReader;
use rust_htslib::bam::{self, record::Cigar};
use std::fs::File;
// use bio::stats::{LogProb, PHREDProb};
use coverm::bam_generator::*;
use estimation::contig_variants::*;
use estimation::variant_matrix::*;
use model::variants::*;
use rayon::prelude::*;
use utils::*;

use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::FlagFilter;
use std::str;

/// Process all reads in a BAM file
#[allow(unused)]
pub fn process_bam<'b, R: IndexedNamedBamReader>(
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
            let quals = record.qual();
            let mut read_cursor: usize = 0;
            let read_len = record.seq().len();
            for cig in record.cigar().iter() {
                match cig {
                    Cigar::Match(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                        // if M, X, or = increment start and decrement end index
                        ups_and_downs[cursor] += 1;
                        let final_pos = cursor + cig.len() as usize;
                        // For checking against mnv
                        let mut potential_mnv = false;
                        let mut mnv_pos = 0;
                        let mut mnv = vec![];
                        let mut mnv_cursor = 0;
                        let mut mnv_qual = 0.;
                        for qpos in read_cursor..(read_cursor + cig.len() as usize) {
                            // See if read is match MNV
                            let qual_pos = record.qual()[qpos] as f64;
                            if potential_mnv && (mnv_pos < mnv.len()) {
                                let read_char = record.seq()[qpos];
                                debug!("MNV searching {} {:?} {}", &mnv_pos, &mnv, &read_char);
                                if mnv[mnv_pos] == read_char {
                                    mnv_pos += 1;
                                    mnv_qual += qual_pos;
                                    debug!("pos {} length {}", &mnv_pos, &mnv.len());
                                    if mnv_pos == mnv.len() {
                                        match variant_matrix.variants(tid, mnv_cursor) {
                                            Some(current_variants) => {
                                                current_variants.iter_mut().for_each(
                                                    |(variant, base)| match variant {
                                                        Variant::MNV(alt) => {
                                                            debug!(
                                                                "alt {:?} found {:?}",
                                                                &alt, &mnv
                                                            );
                                                            if *alt == mnv {
                                                                base.truedepth[sample_idx] += 1;
                                                                base.quals[sample_idx] += qual_pos;
                                                            }
                                                        }
                                                        _ => {
                                                            debug!(
                                                                "Looping through non-MNV variants"
                                                            );
                                                        }
                                                    },
                                                );
                                                mnv = vec![];
                                                mnv_pos = 0;
                                                potential_mnv = false;
                                                mnv_qual = 0.;
                                            }
                                            None => {
                                                debug!("No associated MNV found for {:?}", &mnv);
                                                mnv = vec![];
                                                mnv_pos = 0;
                                                potential_mnv = false;
                                                mnv_qual = 0.;
                                            }
                                        };
                                        mnv = vec![];
                                        mnv_pos = 0;
                                        potential_mnv = false;
                                        mnv_qual = 0.;
                                    }
                                } else {
                                    debug!("Read did not contain correct MNV");
                                    mnv = vec![];
                                    mnv_pos = 0;
                                    potential_mnv = false;
                                    mnv_qual = 0.;
                                }
                            }
                            match variant_matrix.variants(tid, cursor as i64) {
                                Some(current_variants) => {
                                    let read_char = record.seq()[qpos];
                                    current_variants.iter_mut().for_each(|(variant, base)| {
                                        match variant {
                                            Variant::SNV(alt) => {
                                                if *alt == read_char {
                                                    base.truedepth[sample_idx] += 1;
                                                    base.quals[sample_idx] += qual_pos;
                                                }
                                            },
                                            // We need to check every position of the MNV
                                            Variant::MNV(alt) => {
                                                if !potential_mnv {
                                                    if alt[mnv_pos] == read_char {
                                                        mnv = alt.clone();
                                                        debug!("Potential MNV  pos {} var {:?} read {} ref {:?}", &mnv_pos, &mnv, &read_char, &base.refr);

                                                        mnv_pos += 1;
                                                        potential_mnv = true;
                                                        mnv_cursor = cursor as i64;
                                                        mnv_qual = 0.;

                                                        // Then it is automatically assigned
                                                        if mnv_pos == mnv.len() {
                                                            base.truedepth[sample_idx] += 1;
                                                            base.quals[sample_idx] += qual_pos;
                                                            mnv = vec!();
                                                            mnv_pos = 0;
                                                            potential_mnv = false;
                                                            mnv_qual = 0.;
                                                        }
                                                    }
                                                }
                                            },
                                            Variant::None => {
                                                if base.refr[0] == read_char {
                                                    base.truedepth[sample_idx] += 1;
                                                    base.referencedepth[sample_idx] += 1;
                                                    base.quals[sample_idx] += qual_pos;

                                                }
                                            },
                                            _ => {}
                                        }
                                    });
                                }
                                _ => {}
                            }
                            cursor += 1;
                        }

                        // debug!("CIGAR ended, resetting MNV");
                        mnv = vec![];
                        mnv_pos = 0;
                        potential_mnv = false;
                        mnv_qual = 0.;

                        if final_pos < ups_and_downs.len() {
                            // True unless the read hits the contig end.
                            ups_and_downs[final_pos] -= 1;
                        }
                        read_cursor += cig.len() as usize;
                    }
                    Cigar::Del(del) => {
                        match variant_matrix.variants(tid, cursor as i64) {
                            Some(current_variants) => {
                                current_variants.par_iter_mut().for_each(|(variant, base)| {
                                    match variant {
                                        // We need to check every position of the MNV
                                        Variant::Deletion(alt) => {
                                            if alt == del {
                                                base.truedepth[sample_idx] += 1;
                                                base.quals[sample_idx] +=
                                                    record.qual()[read_cursor] as f64
                                            }
                                        }
                                        _ => {}
                                    }
                                });
                            }
                            _ => {}
                        }

                        cursor += cig.len() as usize;
                    }
                    Cigar::RefSkip(del) => {
                        match variant_matrix.variants(tid, cursor as i64) {
                            Some(current_variants) => {
                                current_variants.par_iter_mut().for_each(|(variant, base)| {
                                    match variant {
                                        // We need to check every position of the MNV
                                        Variant::Deletion(alt) => {
                                            if alt == del {
                                                base.truedepth[sample_idx] += 1;
                                                base.quals[sample_idx] +=
                                                    record.qual()[read_cursor] as f64
                                            }
                                        }
                                        _ => {}
                                    }
                                });
                            }
                            _ => {}
                        }
                        cursor += cig.len() as usize;
                    }
                    Cigar::Ins(ins) => {
                        let insertion =
                            &record.seq().as_bytes()[read_cursor..read_cursor + cig.len() as usize];
                        match variant_matrix.variants(tid, cursor as i64) {
                            Some(current_variants) => {
                                current_variants.par_iter_mut().for_each(|(variant, base)| {
                                    match variant {
                                        // We need to check every position of the MNV
                                        Variant::Insertion(alt) => {
                                            if String::from_utf8(alt.to_vec())
                                                .expect("Unable to convert to string")
                                                .contains(
                                                    &String::from_utf8(insertion.to_vec())
                                                        .expect("Unable to convert to string"),
                                                )
                                            {
                                                base.truedepth[sample_idx] += 1;
                                                let qual_sum = quals
                                                    [read_cursor..read_cursor + cig.len() as usize]
                                                    .par_iter()
                                                    .map(|q| *q as f64)
                                                    .sum::<f64>()
                                                    as f64;

                                                base.quals[sample_idx] += qual_sum;
                                            }
                                        }
                                        _ => {}
                                    }
                                });
                            }
                            _ => {}
                        }
                        read_cursor += cig.len() as usize;
                        total_indels_in_current_contig += cig.len() as u64;
                    }
                    Cigar::SoftClip(_) => {
                        // soft clipped portions of reads can be included as structural variants
                        // not sure if this correct protocol or not
                        let insertion =
                            &record.seq().as_bytes()[read_cursor..read_cursor + cig.len() as usize];
                        match variant_matrix.variants(tid, cursor as i64) {
                            Some(current_variants) => {
                                current_variants.par_iter_mut().for_each(|(variant, base)| {
                                    match variant {
                                        // We need to check every position of the MNV
                                        Variant::Insertion(alt) | Variant::Inversion(alt) => {
                                            if String::from_utf8(alt.to_vec())
                                                .expect("Unable to convert to string")
                                                .contains(
                                                    &String::from_utf8(insertion.to_vec())
                                                        .expect("Unable to convert to string"),
                                                )
                                            {
                                                base.truedepth[sample_idx] += 1;
                                                let qual_sum = quals
                                                    [read_cursor..read_cursor + cig.len() as usize]
                                                    .par_iter()
                                                    .map(|q| *q as f64)
                                                    .sum::<f64>()
                                                    as f64;

                                                base.quals[sample_idx] += qual_sum;
                                            }
                                        }
                                        _ => {}
                                    }
                                });
                            }
                            _ => {}
                        }
                        read_cursor += cig.len() as usize;
                    }
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
        tid as i32,
        ups_and_downs,
        coverage_estimators,
        min,
        max,
        total_indels_in_current_contig as usize,
        contig_end_exclusion,
        target_len as usize,
        contig_name,
        variant_matrix,
        sample_idx,
        method,
        total_mismatches,
        coverage_fold,
        num_mapped_reads_in_current_contig,
        sample_count,
        output_prefix,
        &stoit_name,
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

    // if bam_generated.num_detected_primary_alignments() == 0 {
    //     warn!(
    //         "No primary alignments were observed for sample {} for {} \
    //        - perhaps something went wrong in the mapping?",
    //         stoit_name,
    //         &reference,
    //     );
    // }
    bam_generated.finish();
}

#[allow(unused)]
pub fn process_previous_contigs_var(
    mode: &str,
    last_tid: i32,
    ups_and_downs: Vec<i32>,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    min: f32,
    max: f32,
    total_indels_in_current_contig: usize,
    contig_end_exclusion: u64,
    contig_len: usize,
    contig_name: Vec<u8>,
    variant_matrix: &mut VariantMatrix,
    sample_idx: usize,
    method: &str,
    total_mismatches: u64,
    coverage_fold: f32,
    num_mapped_reads_in_current_contig: u64,
    sample_count: usize,
    output_prefix: &str,
    stoit_name: &str,
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

        let mut variant_struct =
            VariantStats::new_contig_stats(min as f64, max as f64, contig_end_exclusion);

        // adds contig info to variant struct
        variant_struct.add_contig(
            variant_matrix.variants_of_contig(last_tid),
            last_tid.clone(),
            total_indels_in_current_contig,
            contig_name.clone(),
            contig_len,
            sample_idx,
            coverages,
            ups_and_downs,
        );

        match mode {
            "bin" => {
                // Add samples contig information to main struct
                debug!("Adding in new info for contig...");
                variant_matrix.add_contig(variant_struct, sample_count, sample_idx);
            }
            _ => {
                panic!("unknown mode {}", mode);
            }
        }
    }
}
