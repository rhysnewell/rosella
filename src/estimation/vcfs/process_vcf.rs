use bio::io::fasta::IndexedReader;
use bird_tool_utils::command;
use rust_htslib::errors::Error;
use rust_htslib::{bcf, bcf::Read};
use std;
use std::collections::{HashMap, HashSet};
use std::fs::File;

use coverm::bam_generator::*;
use coverm::FlagFilter;
use estimation::variant_matrix::*;
use external_command_checker;
use model::variants::*;
use utils::*;

use crate::*;
use scoped_threadpool::Pool;
use statrs::statistics::Variance;
use std::io::Write;
use std::path::Path;
use std::str;
use std::sync::{Arc, Mutex};
use tempdir::TempDir;
use tempfile::Builder;

#[allow(unused)]
pub fn process_vcf<'b, R: IndexedNamedBamReader + Send, G: NamedBamReaderGenerator<R> + Send>(
    bam_generator: G,
    split_threads: usize,
    sample_idx: usize,
    mut sample_count: usize,
    variant_matrix: &'b mut VariantMatrix,
    readtype: ReadType,
    m: &'b clap::ArgMatches,
    reference_file: &mut IndexedReader<File>,
    reference: &str,
    mut short_sample_count: usize,
    flag_filters: &'b FlagFilter,
    tid: u32,
) {
    let mut bam_generated = bam_generator.start();
    let mut stoit_name = bam_generated.name().to_string().replace("/", ".");

    bam_generated.set_threads(split_threads);

    let header = bam_generated.header().clone(); // bam header
    let target_len = header.target_len(tid).unwrap();
    let target_name = header.target_names()[tid as usize];
    let target_name = String::from_utf8(target_name.to_vec()).unwrap();

    let bam_path = bam_generated.path().to_string();

    // minimum PHRED base quality
    let bq = m
        .value_of("base-quality-threshold")
        .unwrap()
        .parse::<f64>()
        .unwrap();
    // Minimum MAPQ value
    let mapq_thresh = std::cmp::max(1, m.value_of("mapq-threshold").unwrap().parse().unwrap());
    // Minimum count for a SNP to be considered
    let min_variant_depth: i32 = m.value_of("min-variant-depth").unwrap().parse().unwrap();
    let min_variant_quality = m
        .value_of("min-variant-quality")
        .unwrap()
        .parse::<f64>()
        .unwrap();

    let mut total_records = 0;
    variant_matrix.add_sample_name(stoit_name.to_string(), sample_idx);
    let mut contig_stats: Vec<f64> = Vec::new();
    let valid_chars: HashSet<u8> = vec![
        "A".as_bytes()[0],
        "T".as_bytes()[0],
        "C".as_bytes()[0],
        "G".as_bytes()[0],
    ]
    .into_iter()
    .collect::<HashSet<u8>>();
    // for each genomic position, only has hashmap when variants are present. Includes read ids
    match readtype {
        ReadType::Short | ReadType::Long => {
            // target_names
            //     .iter()
            //     .enumerate()
            //     .for_each(|(tid, contig_name)| {

            {
                // use pileups to call SNPs for low quality variants
                // That are usually skipped by GATK
                let mut coverage = Vec::with_capacity(target_len as usize);
                let mut depth_sum = 0;
                variant_matrix.add_info(tid as usize, target_name.as_bytes().to_vec(), target_len);
                {
                    bam_generated.fetch((tid as u32));
                    match bam_generated.pileup() {
                        Some(pileups) => {
                            let mut ref_seq = Vec::new();
                            // Update all contig information
                            fetch_contig_from_reference(
                                reference_file,
                                &target_name.as_bytes().to_vec(),
                            );
                            read_sequence_to_vec(
                                &mut ref_seq,
                                reference_file,
                                &target_name.as_bytes().to_vec(),
                            );

                            for p in pileups {
                                // if {
                                let pileup = p.unwrap();
                                let tid = pileup.tid() as i32;
                                let pos = pileup.pos() as usize;
                                let depth = pileup.depth();
                                coverage.push(depth as f64);
                                depth_sum += depth;
                                let refr_base = ref_seq[pos];

                                if refr_base != "N".as_bytes()[0] {
                                    // let mut base = Base::new(tid, pos, );
                                    // info!("Base dict {:?} {}", &pos_dict, pos);

                                    let mut base_dict = HashMap::new();

                                    let mut refr_depth = 0;
                                    let mut refr_qual = 0.;

                                    for alignment in pileup.alignments() {
                                        let record = alignment.record();
                                        if (!flag_filters.include_supplementary
                                            && record.is_supplementary()
                                            && readtype != ReadType::Long)
                                            || (!flag_filters.include_secondary
                                                && record.is_secondary())
                                            || (!flag_filters.include_improper_pairs
                                                && !record.is_proper_pair()
                                                && readtype != ReadType::Long)
                                        {
                                            continue;
                                        } else if !flag_filters.include_secondary
                                            && record.is_secondary()
                                            && readtype == ReadType::Long
                                        {
                                            continue;
                                        }
                                        if record.mapq() >= mapq_thresh {
                                            if !alignment.is_del() && !alignment.is_refskip() {
                                                // query position in read
                                                let qpos = alignment.qpos().unwrap();
                                                let record_qual = record.qual()[qpos] as f64;
                                                if record_qual >= bq {
                                                    let read_base = alignment.record().seq()[qpos];
                                                    if read_base != refr_base
                                                        && valid_chars.contains(&read_base)
                                                    {
                                                        let mut base = base_dict
                                                            .entry(read_base)
                                                            .or_insert(Base::new(
                                                                tid as u32,
                                                                pos as i64,
                                                                sample_count,
                                                                vec![refr_base],
                                                            ));

                                                        base.depth[sample_idx] += 1;
                                                        base.quals[sample_idx] += record_qual;
                                                        if base.variant == Variant::None {
                                                            base.variant = Variant::SNV(read_base);
                                                        }
                                                    }
                                                } else {
                                                    refr_depth += 1;
                                                    refr_qual += record_qual;
                                                }
                                            }
                                        }
                                    }

                                    // Collect refr base information
                                    {
                                        let mut base =
                                            base_dict.entry(refr_base).or_insert(Base::new(
                                                tid as u32,
                                                pos as i64,
                                                sample_count,
                                                vec![refr_base],
                                            ));

                                        base.depth[sample_idx] = refr_depth;
                                        base.quals[sample_idx] = refr_qual;
                                    }

                                    // If more than one variant at location (including reference)
                                    // Collect the variants
                                    if base_dict.keys().len() > 1 {
                                        let mut variant_found = false;
                                        let mut refr_base = 0;
                                        for (var_char, base) in base_dict.iter() {
                                            if base.depth[sample_idx] >= min_variant_depth &&
                                                // base.quals[sample_idx] >= min_variant_quality &&
                                                base.variant != Variant::None
                                            {
                                                total_records += 1;
                                                variant_found = true;
                                                refr_base = base.refr[0];
                                                variant_matrix.add_variant_to_matrix(
                                                    sample_idx,
                                                    base,
                                                    tid as usize,
                                                );
                                            }
                                        }
                                        // Add reference
                                        match base_dict.get(&refr_base) {
                                            Some(base) => {
                                                variant_matrix.add_variant_to_matrix(
                                                    sample_idx,
                                                    base,
                                                    tid as usize,
                                                );
                                            }
                                            None => {}
                                        }
                                    }
                                }
                            }
                        }
                        None => println!("no bam for pileups"),
                    };
                    // bam_generated.set_threads(1);
                };
                // Get coverage mean and standard dev
                let contig_cov = depth_sum as f64 / target_len as f64;
                let std_dev = 10. * coverage.std_dev();
                contig_stats = vec![contig_cov, std_dev];
                // });
            }
        }
        _ => {}
    }

    bam_generated.finish();
    let mut variant_matrix_sync = Arc::new(Mutex::new(variant_matrix.clone()));
    let freebayes_threads = std::cmp::max(split_threads, 1);

    if readtype == ReadType::Long
        || readtype == ReadType::Assembly
        || (readtype == ReadType::Short && m.is_present("freebayes"))
    {
        //placeholder check
        // Get VCF file from BAM using freebayes of SVIM
        let mut vcf_reader = get_vcf(
            &stoit_name,
            &m,
            freebayes_threads,
            &readtype,
            target_len,
            &reference,
            &bam_path,
            &target_name,
            tid,
        );

        match vcf_reader {
            Ok(ref mut reader) => {
                let min_qual = m.value_of("min-variant-quality").unwrap().parse().unwrap();

                let mut pool = Pool::new(freebayes_threads as u32);

                pool.scoped(|scope| {
                    for vcf_record in reader.records().into_iter() {
                        scope.execute(|| {
                            let mut vcf_record = vcf_record.unwrap();
                            let variant_rid = vcf_record.rid().unwrap();

                            let base_option = Base::from_vcf_record(
                                &mut vcf_record,
                                sample_count,
                                sample_idx,
                                &readtype,
                                min_qual,
                            );
                            match base_option {
                                Some(bases) => {
                                    for base in bases {
                                        let mut variant_matrix_sync =
                                            variant_matrix_sync.lock().unwrap();
                                        variant_matrix_sync.add_variant_to_matrix(
                                            sample_idx,
                                            &base,
                                            variant_rid as usize,
                                        );
                                    }
                                }
                                None => {}
                            }
                        });
                    }
                });
            }
            Err(_) => {
                debug!("No VCF records found for sample {}", &stoit_name);
            }
        }
    }

    let mut variant_matrix_sync = variant_matrix_sync.lock().unwrap();
    if readtype == ReadType::Short || readtype == ReadType::Long {
        variant_matrix_sync.remove_variants(tid as i32, sample_idx, contig_stats);
    }
    //     }
    // });
    *variant_matrix = variant_matrix_sync.clone();
}

/// Get or generate vcf file
#[allow(unused)]
pub fn get_vcf(
    stoit_name: &str,
    m: &clap::ArgMatches,
    threads: usize,
    readtype: &ReadType,
    target_length: u64,
    reference: &str,
    bam_path: &str,
    target_name: &String,
    tid: u32,
) -> std::result::Result<bcf::Reader, Error> {
    return generate_vcf(
        bam_path,
        m,
        threads,
        readtype,
        target_length,
        reference,
        target_name,
        tid,
    );
}

/// Makes direct call to freebayes or SVIM
#[allow(unused)]
pub fn generate_vcf(
    bam_path: &str,
    m: &clap::ArgMatches,
    threads: usize,
    readtype: &ReadType,
    target_length: u64,
    reference: &str,
    target_name: &String,
    tid: u32,
) -> std::result::Result<bcf::Reader, Error> {
    // setup temp directory
    let tmp_dir = TempDir::new("lorikeet_fifo").expect("Unable to create temporary directory");

    let tmp_bam_path1 = Builder::new()
        .prefix(&(tmp_dir.path().to_str().unwrap().to_string() + "/"))
        .suffix(".bam")
        .tempfile()
        .unwrap();

    let mapq_thresh = std::cmp::max(1, m.value_of("mapq-threshold").unwrap().parse().unwrap());

    match readtype {
        &ReadType::Short => {
            external_command_checker::check_for_freebayes();
            external_command_checker::check_for_freebayes_parallel();
            external_command_checker::check_for_samtools();
            external_command_checker::check_for_vt();

            let region_size = 100000;

            let index_path = format!("{}.fai", reference);

            let vcf_path = &(tmp_dir.path().to_str().unwrap().to_string() + "/output.vcf");
            let vcf_path_prenormalization =
                &(tmp_dir.path().to_str().unwrap().to_string() + "/output_prenormalization.vcf");

            let mut region_tmp_file = Builder::new()
                .prefix(&(tmp_dir.path().to_str().unwrap().to_string() + "/"))
                .suffix(".txt")
                .tempfile()
                .unwrap();

            let mut total_region_covered = 0;
            let total_region_chunks = target_length / region_size
                + if target_length % region_size != 0 {
                    1
                } else {
                    0
                };
            for idx in (0..total_region_chunks).into_iter() {
                total_region_covered += region_size;
                if total_region_covered > target_length {
                    writeln!(
                        region_tmp_file,
                        "{}:{}-{}",
                        &target_name,
                        total_region_covered - region_size,
                        target_length,
                    )
                    .unwrap();
                } else {
                    writeln!(
                        region_tmp_file,
                        "{}:{}-{}",
                        &target_name,
                        total_region_covered - region_size,
                        total_region_covered,
                    )
                    .unwrap();
                }
            }

            // Variant calling pipeline adapted from Snippy but without all of the rewriting of BAM files
            let vcf_cmd_string = format!(
                "set -e -o pipefail;  \
                    ulimit -s {} && freebayes-parallel {:?} {} -f {} -C {} -q {} \
                    --min-repeat-entropy {} -p {} --strict-vcf -m {} {} > {}",
                m.value_of("ulimit").unwrap(),
                region_tmp_file.path(),
                threads,
                &reference,
                m.value_of("min-variant-depth").unwrap(),
                m.value_of("base-quality-threshold").unwrap(),
                m.value_of("min-repeat-entropy").unwrap(),
                m.value_of("ploidy").unwrap(),
                mapq_thresh,
                &bam_path,
                &vcf_path_prenormalization,
            );
            let vt_cmd_string = format!(
                "vt normalize {} -n -r {} > {}",
                &vcf_path_prenormalization,
                &reference,
                &vcf_path,
                // &vcf_path,
            );
            debug!("Queuing cmd_string: {}", vcf_cmd_string);
            command::finish_command_safely(
                std::process::Command::new("bash")
                    .arg("-c")
                    .arg(&vcf_cmd_string)
                    .stderr(std::process::Stdio::piped())
                    // .stdout(std::process::Stdio::piped())
                    .spawn()
                    .expect("Unable to execute bash"),
                "freebayes",
            );
            debug!("Queuing cmd_string: {}", vt_cmd_string);
            command::finish_command_safely(
                std::process::Command::new("bash")
                    .arg("-c")
                    .arg(&vt_cmd_string)
                    .stderr(std::process::Stdio::piped())
                    .stdout(std::process::Stdio::piped())
                    .spawn()
                    .expect("Unable to execute bash"),
                "vt",
            );
            debug!("VCF Path {:?}", vcf_path);
            let vcf_reader = bcf::Reader::from_path(&Path::new(vcf_path));

            tmp_dir.close().expect("Failed to close temp directory");
            return vcf_reader;
        }
        &ReadType::Long => {
            external_command_checker::check_for_svim();
            let svim_path = tmp_dir.path().to_str().unwrap().to_string();

            let cmd_string = format!(
                "set -e -o pipefail; samtools view -bh {} {} > {} && \
        svim alignment --read_names --skip_genotyping \
        --tandem_duplications_as_insertions --interspersed_duplications_as_insertions \
        --min_mapq {} --sequence_alleles {} {} {}",
                bam_path,
                &target_name,
                tmp_bam_path1.path().to_str().unwrap().to_string(),
                mapq_thresh,
                &svim_path,
                tmp_bam_path1.path().to_str().unwrap().to_string(),
                &reference
            );
            debug!("Queuing cmd_string: {}", cmd_string);
            command::finish_command_safely(
                std::process::Command::new("bash")
                    .arg("-c")
                    .arg(&cmd_string)
                    .stderr(std::process::Stdio::piped())
                    //                .stdout(std::process::Stdio::null())
                    .spawn()
                    .expect("Unable to execute bash"),
                "svim",
            );
            let vcf_path = &(svim_path + "/variants.vcf");
            debug!("VCF Path {:?}", vcf_path);
            let vcf_reader = bcf::Reader::from_path(&Path::new(vcf_path));

            tmp_dir.close().expect("Failed to close temp directory");
            return vcf_reader;
        }
        &ReadType::Assembly => {
            external_command_checker::check_for_svim_asm();
            let svim_path = tmp_dir.path().to_str().unwrap().to_string();

            let cmd_string = format!(
                "set -e -o pipefail; samtools view -bh {} {} > {} && \
        svim-asm haploid --query_names \
        --tandem_duplications_as_insertions --interspersed_duplications_as_insertions \
        --min_mapq {} {} {} {}",
                bam_path,
                &target_name,
                tmp_bam_path1.path().to_str().unwrap().to_string(),
                mapq_thresh,
                &svim_path,
                tmp_bam_path1.path().to_str().unwrap().to_string(),
                &reference
            );
            debug!("Queuing cmd_string: {}", cmd_string);
            command::finish_command_safely(
                std::process::Command::new("bash")
                    .arg("-c")
                    .arg(&cmd_string)
                    .stderr(std::process::Stdio::piped())
                    //                .stdout(std::process::Stdio::null())
                    .spawn()
                    .expect("Unable to execute bash"),
                "svim",
            );
            let vcf_path = &(svim_path + "/variants.vcf");
            debug!("VCF Path {:?}", vcf_path);
            let vcf_reader = bcf::Reader::from_path(&Path::new(vcf_path));

            tmp_dir.close().expect("Failed to close temp directory");
            return vcf_reader;
        }
    }
}
