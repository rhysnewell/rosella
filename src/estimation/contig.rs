use std;
use std::collections::HashMap;

use bird_tool_utils::command;
use coverm::bam_generator::*;
use estimation::bams::{index_bams::*, process_bam::*};
use estimation::variant_matrix::*;
use estimation::vcfs::process_vcf::*;
use external_command_checker;
use utils::*;

use crate::*;
use bio::io::gff;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::FlagFilter;
use glob::glob;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rayon::prelude::*;
use scoped_threadpool::Pool;
use std::path::Path;
use std::process::Stdio;
use std::str;
use std::sync::{Arc, Mutex};
use tempdir::TempDir;
use tempfile::NamedTempFile;

#[derive(Clone, Debug)]
struct Elem {
    key: String,
    index: usize,
    progress_bar: ProgressBar,
}

#[allow(unused)]
pub fn pileup_variants<
    'a,
    R: NamedBamReader,
    S: NamedBamReaderGenerator<R>,
    T: NamedBamReader,
    U: NamedBamReaderGenerator<T>,
    V: NamedBamReader,
    W: NamedBamReaderGenerator<V>,
>(
    m: &clap::ArgMatches,
    bam_readers: Vec<S>,
    longreads: Option<Vec<U>>,
    assembly_readers: Option<Vec<W>>,
    mode: &str,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    flag_filters: FlagFilter,
    mapq_threshold: u8,
    min_var_depth: usize,
    min: f32,
    max: f32,
    contig_end_exclusion: u64,
    output_prefix: &str,
    n_threads: usize,
    method: &str,
    coverage_fold: f32,
    include_indels: bool,
    include_soft_clipping: bool,
    is_long_read: bool,
    genomes_and_contigs: GenomesAndContigs,
    tmp_bam_file_cache: Option<TempDir>,
    concatenated_genomes: Option<NamedTempFile>,
) {
    // TODO: Split up analysis per contig for speed purposes
    let reference = m.value_of("reference").unwrap();
    let mut indexed_reference = generate_faidx(reference);
    // All different counts of samples I need. Changes depends on when using concatenated genomes or not
    let mut short_sample_count = bam_readers.len();
    let mut long_sample_count = 0;
    let mut assembly_sample_count = 0;

    let mut ani = 0.;

    let longreads = match longreads {
        Some(vec) => {
            long_sample_count += vec.len();
            vec
        }
        None => vec![],
    };

    let assembly = match assembly_readers {
        Some(vec) => {
            assembly_sample_count += vec.len();
            vec
        }
        None => vec![],
    };

    let alpha: f64 = m.value_of("fdr-threshold").unwrap().parse().unwrap();

    // Finish each BAM source
    if m.is_present("longreads") || m.is_present("longread-bam-files") {
        info!("Processing long reads...");
        finish_bams(longreads, n_threads);
    }

    if m.is_present("assembly") || m.is_present("assembly-bam-files") {
        info!("Processing assembly alignments...");
        finish_bams(assembly, n_threads);
    }
    // if !m.is_present("bam-files") {
    info!("Processing short reads...");
    finish_bams(bam_readers, n_threads);
    // }

    // Put reference index in the variant map and initialize matrix
    let mut progress_bars = vec![
        Elem {
            key: "Genomes complete".to_string(),
            index: 1,
            progress_bar: ProgressBar::new(0),
        };
        2
    ];


    progress_bars[0] = Elem {
        key: "Operations remaining".to_string(),
        index: 0,
        progress_bar: ProgressBar::new(
            (short_sample_count + long_sample_count + assembly_sample_count) as u64
                * 2
                + 1,
        ),
    };

    debug!(
        "{} Longread BAM files, {} Shortread BAM files and {} assembly alignment BAMs {} Total BAMs",
        long_sample_count,
        short_sample_count,
        assembly_sample_count,
        (short_sample_count + long_sample_count + assembly_sample_count)
    );

    let parallel_genomes = m.value_of("parallel-genomes").unwrap().parse().unwrap();
    let mut pool = Pool::new(parallel_genomes);
    let n_threads = std::cmp::max(n_threads / parallel_genomes as usize, 2);
    // Set up multi progress bars
    let multi = Arc::new(MultiProgress::new());
    let sty_eta = ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg} ETA: [{eta}]");

    let sty_aux = ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {spinner:.green} {msg} {pos:>4}/{len:4}");
    progress_bars
        .par_iter()
        .for_each(|pb| pb.progress_bar.set_style(sty_aux.clone()));
    progress_bars[0].progress_bar.set_style(sty_eta.clone());

    // let pb_main = multi.add(ProgressBar::new(reference_map.keys().len() as u64));
    // pb_main.set_style(sty_eta.clone());

    let tree: Arc<Mutex<Vec<&Elem>>> =
        Arc::new(Mutex::new(Vec::with_capacity(progress_bars.len())));
    {
        let mut tree = tree.lock().unwrap();
        for pb in progress_bars.iter() {
            tree.push(pb)
        }
    }
    // let tree2 = Arc::clone(&tree);

    let multi_inner = Arc::clone(&multi);

    pool.scoped(|scope| {
        {
            // Total steps eta progress bar
            let elem = &progress_bars[0];
            let pb = multi_inner.insert(0, elem.progress_bar.clone());

            pb.enable_steady_tick(500);

            pb.set_message(&format!("{}...", &elem.key,));
        }
        {
            // completed genomes progress bar
            let elem = &progress_bars[1];
            let pb = multi_inner.insert(1, elem.progress_bar.clone());

            pb.enable_steady_tick(500);

            pb.set_message(&format!("{}...", &elem.key,));
        }
            let multi_inner = &multi_inner;
            let tree = &tree;
            let progress_bars = &progress_bars;
            let flag_filters = &flag_filters;
            let tmp_bam_file_cache = match tmp_bam_file_cache.as_ref() {
                Some(cache) => Some(cache.path().to_str().unwrap().to_string()),
                None => None,
            };

            let mut coverage_estimators = coverage_estimators.clone();

            // pb_main.tick();
            //
            // pb_main.inc(1);
            // pb_main.set_message(&format!(
            //     "Staging reference: {}",
            //     &genomes_and_contigs.genomes[ref_idx],
            // ));

            scope.execute(move || {

                {
                    let elem = &progress_bars[2];
                    let pb = multi_inner.insert(2, elem.progress_bar.clone());

                    pb.enable_steady_tick(500);

                    pb.set_message(&format!("{}: Preparing variants...", &elem.key,));
                    // multi.join().unwrap();

                    // tree.lock().unwrap().insert(elem.index, &elem);
                }

                // Read BAMs back in as indexed
                let mut indexed_bam_readers = recover_bams(
                    m,
                    short_sample_count,
                    long_sample_count,
                    assembly_sample_count,
                    &tmp_bam_file_cache,
                );
                let mut per_reference_samples = 0;
                let mut per_reference_short_samples = 0;
                let mut variant_matrix = VariantMatrix::new_matrix(short_sample_count + long_sample_count + assembly_sample_count);

                debug!(
                    "Running SNP calling on {} samples",
                    indexed_bam_readers.len()
                );


                indexed_bam_readers.into_iter().enumerate().for_each(
                    |(sample_idx, bam_generator)| {
                        // Get the appropriate sample index based on how many references we are using
                        let mut bam_generator = generate_indexed_named_bam_readers_from_bam_files(
                            vec![&bam_generator],
                            n_threads as u32,
                        )
                        .into_iter()
                        .next()
                        .unwrap();
                        if sample_idx < short_sample_count {
                            process_vcf(
                                bam_generator,
                                n_threads,
                                sample_idx,
                                per_reference_samples,
                                &mut variant_matrix,
                                ReadType::Short,
                                m,
                                &mut indexed_reference,
                                &reference,
                                per_reference_short_samples,
                                &flag_filters,
                            );
                        } else if (m.is_present("longreads") | m.is_present("longread-bam-files"))
                            && sample_idx >= short_sample_count
                            && sample_idx < (short_sample_count + long_sample_count)
                        {
                            debug!("Running structural variant detection...");
                            // Get the appropriate sample index based on how many references we are using by tracking
                            // changes in references
                            process_vcf(
                                bam_generator,
                                n_threads,
                                sample_idx,
                                per_reference_samples,
                                &mut variant_matrix,
                                ReadType::Long,
                                m,
                                &mut indexed_reference,
                                &reference,
                                per_reference_short_samples,
                                &flag_filters,
                            );
                        } else if (m.is_present("assembly") | m.is_present("assembly_bam_files"))
                            && sample_idx >= (short_sample_count + long_sample_count)
                        {
                            process_vcf(
                                bam_generator,
                                n_threads,
                                sample_idx,
                                per_reference_samples,
                                &mut variant_matrix,
                                ReadType::Assembly,
                                m,
                                &mut indexed_reference,
                                &reference,
                                per_reference_short_samples,
                                &flag_filters,
                            );
                        }
                        {
                            let pb = &tree.lock().unwrap()[2];

                            pb.progress_bar.set_message(&format!(
                                "{}: Variant calling on sample: {}",
                                pb.key,
                                variant_matrix.get_sample_name(sample_idx),
                            ));
                            pb.progress_bar.inc(1);
                        }
                        {
                            let pb = &tree.lock().unwrap()[0];
                            pb.progress_bar.inc(1);
                        }
                    },
                );
                {
                    let pb = &tree.lock().unwrap()[2];
                    pb.progress_bar
                        .set_message(&format!("{}: Initial variant calling complete...", pb.key));
                }
                // // Read BAMs back in as indexed
                let mut indexed_bam_readers = recover_bams(
                    m,
                    short_sample_count,
                    long_sample_count,
                    assembly_sample_count,
                    &tmp_bam_file_cache,
                );

                {
                    let pb = &tree.lock().unwrap()[2];
                    pb.progress_bar.reset();
                    pb.progress_bar.enable_steady_tick(1000);
                    pb.progress_bar
                        .set_message(&format!("{}: Performing guided variant calling...", pb.key));
                }
                // let mut variant_matrix = Mutex::new(variant_matrix);
                if variant_matrix.get_variant_count() > 0 {
                    // if there are variants, perform guided variant calling
                    std::fs::create_dir_all(&output_prefix).unwrap();

                    indexed_bam_readers.into_iter().enumerate().for_each(
                        |(sample_idx, bam_generator)| {
                            let mut bam_generator =
                                generate_indexed_named_bam_readers_from_bam_files(
                                    vec![&bam_generator],
                                    n_threads as u32,
                                )
                                .into_iter()
                                .next()
                                .unwrap();
                            if sample_idx < short_sample_count {
                                process_bam(
                                    bam_generator,
                                    sample_idx,
                                    per_reference_samples,
                                    &mut coverage_estimators,
                                    &mut variant_matrix,
                                    n_threads,
                                    m,
                                    &output_prefix,
                                    coverage_fold,
                                    min_var_depth,
                                    contig_end_exclusion,
                                    min,
                                    max,
                                    mode,
                                    include_soft_clipping,
                                    include_indels,
                                    &flag_filters,
                                    mapq_threshold,
                                    method,
                                    ReadType::Short,
                                    &mut indexed_reference,
                                )
                            } else if sample_idx >= short_sample_count
                                && sample_idx < (short_sample_count + long_sample_count)
                            {
                                process_bam(
                                    bam_generator,
                                    sample_idx,
                                    per_reference_samples,
                                    &mut coverage_estimators,
                                    &mut variant_matrix,
                                    n_threads,
                                    m,
                                    &output_prefix,
                                    coverage_fold,
                                    min_var_depth,
                                    contig_end_exclusion,
                                    min,
                                    max,
                                    mode,
                                    include_soft_clipping,
                                    include_indels,
                                    &flag_filters,
                                    mapq_threshold,
                                    method,
                                    ReadType::Long,
                                    &mut indexed_reference,
                                )
                            } else if sample_idx >= (short_sample_count + long_sample_count) {
                                process_bam(
                                    bam_generator,
                                    sample_idx,
                                    per_reference_samples,
                                    &mut coverage_estimators,
                                    &mut variant_matrix,
                                    n_threads,
                                    m,
                                    &output_prefix,
                                    coverage_fold,
                                    min_var_depth,
                                    contig_end_exclusion,
                                    min,
                                    max,
                                    mode,
                                    include_soft_clipping,
                                    include_indels,
                                    &flag_filters,
                                    mapq_threshold,
                                    method,
                                    ReadType::Assembly,
                                    &mut indexed_reference,
                                )
                            }

                            {
                                let pb = &tree.lock().unwrap()[2];
                                pb.progress_bar.set_message(&format!(
                                    "{}: Guided variant calling on sample: {}",
                                    pb.key,
                                    variant_matrix.get_sample_name(sample_idx),
                                ));
                                pb.progress_bar.inc(1);
                            }
                            {
                                let pb = &tree.lock().unwrap()[0];
                                pb.progress_bar.inc(1);
                            }
                        },
                    );
                } else {
                    let pb = &tree.lock().unwrap()[0];
                    pb.progress_bar.inc(indexed_bam_readers.len() as u64);
                    pb.progress_bar.reset_eta();
                }
                {
                    let pb = &tree.lock().unwrap()[2];
                    pb.progress_bar
                        .set_message(&format!("{}: Guided variant calling complete...", pb.key));
                }

                // Collects info about variants across samples to check whether they are genuine or not
                // using FDR

                // TODO: Make sure that this is fixed. It seems to work appropriately now
                {
                    let pb = &tree.lock().unwrap()[2];
                    pb.progress_bar
                        .set_message(&format!("{}: Setting FDR threshold...", pb.key));
                }
                variant_matrix
                    .remove_false_discoveries(alpha, reference);

                if mode == "bin" {

                    let anchor_size: usize = m.value_of("n-neighbors").unwrap().parse().unwrap();
                    let n_components: usize = m.value_of("n-components").unwrap().parse().unwrap();

                    // Calculate the geometric mean values and CLR for each variant, reference specific
                    {
                        let pb = &tree.lock().unwrap()[2];
                        pb.progress_bar.set_message(&format!(
                            "{}: Generating variant distances...",
                            &reference,
                        ));
                    }
                    variant_matrix.generate_distances();

                    // Generate initial read linked clusters
                    // Cluster each variant using phi-D and fuzzy DBSCAN, reference specific
                    {
                        let pb = &tree.lock().unwrap()[2];
                        pb.progress_bar
                            .set_message(&format!("{}: Running UMAP and HDBSCAN...", &reference,));
                    }

                    // Get sample distances
                    {
                        let pb = &tree.lock().unwrap()[2];
                        pb.progress_bar.set_message(&format!(
                            "{}: Generating adjacency matrix...",
                            &reference,
                        ));
                    }
                    variant_matrix.calculate_sample_distances(
                        &output_prefix,
                    );

                    // If flagged, then create plots using CMplot
                    let window_size = m.value_of("window-size").unwrap().parse().unwrap();
                    variant_matrix.print_variant_stats(
                        window_size,
                        &output_prefix,
                    );

                    // Write variants in VCF format, reference specific
                    {
                        let pb = &tree.lock().unwrap()[2];
                        pb.progress_bar
                            .set_message(&format!("{}: Generating VCF file...", &reference,));
                    }
                    variant_matrix.write_vcf(&output_prefix,);
                };
                {
                    let pb = &tree.lock().unwrap()[2];
                    pb.progress_bar
                        .set_message(&format!("{}: All steps completed {}", &reference, "✔",));
                    pb.progress_bar.finish_and_clear();
                }
                {
                    let pb = &tree.lock().unwrap()[1];
                    pb.progress_bar.inc(1);
                    let pos = pb.progress_bar.position();
                    let len = pb.progress_bar.length();
                    if pos >= len {
                        pb.progress_bar
                            .finish_with_message(&format!("All genomes analyzed {}", "✔",));
                    }
                }
                {
                    let pb = &tree.lock().unwrap()[0];
                    pb.progress_bar.inc(1);
                    let pos = pb.progress_bar.position();
                    let len = pb.progress_bar.length();
                    if pos >= len {
                        pb.progress_bar
                            .finish_with_message(&format!("All steps completed {}", "✔",));
                    }
                }
            });
        // pb_main.finish_with_message("All genomes staged...");
        multi.join().unwrap();
    });

    info!("Analysis finished!");
    // multi.join_and_clear().unwrap();
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
    //
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
    //            reads_mapped_vec = variant_variants(
    //                bam_readers,
    //                &mut coverage_taker,
    //                coverage_estimators,
    //                print_zero_coverage_contigs,
    //                flag_filters,
    //                false,
    //            );
    //        }
    ////        assert_eq!(expected, str::from_utf8(stream.get_ref()).unwrap());
    //    }
}
