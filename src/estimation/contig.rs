use coverm::bam_generator::*;
use estimation::bams::{contig_coverage::*, index_bams::*, process_bam::*};
use estimation::variant_matrix::*;
use std;
use utils::*;

use crate::*;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::FlagFilter;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use scoped_threadpool::Pool;
use std::collections::HashMap;
use std::path::Path;
use std::str;
use std::sync::Arc;
use tempdir::TempDir;

#[derive(Clone, Debug)]
struct Elem<'a> {
    key: &'a str,
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
    bam_readers: Option<Vec<S>>,
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
    tmp_bam_file_cache: Option<TempDir>,
) {
    // TODO: Split up analysis per contig for speed purposes
    let reference = m.value_of("reference").unwrap();
    generate_faidx(reference); // Check if faidx is present
                               // get number of contigs from faidx
    let n_contigs = bio::io::fasta::Index::from_file(&format!("{}.fai", &reference))
        .unwrap()
        .sequences()
        .len();
    // All different counts of samples I need. Changes depends on when using concatenated genomes or not
    let mut short_sample_count = 0;
    let mut long_sample_count = 0;
    let mut assembly_sample_count = 0;

    let mut ani = 0.;

    let bam_readers = match bam_readers {
        Some(vec) => {
            short_sample_count += vec.len();
            vec
        }
        None => vec![],
    };

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

    // Variant matrix that will collect results from the per contig variant matrices
    // Does not include assembly alignments as samples
    let mut main_variant_matrix = VariantMatrix::new_matrix(
        short_sample_count + long_sample_count,
        Some(bio::io::fasta::Index::from_file(&format!("{}.fai", &reference)).unwrap()),
    );

    let alpha: f64 = m.value_of("fdr-threshold").unwrap().parse().unwrap();

    // }

    // Put reference index in the variant map and initialize matrix
    let mut progress_bars = vec![
        Elem {
            key: "Contigs analyzed",
            index: 1,
            progress_bar: ProgressBar::new(n_contigs as u64),
        };
        1
    ];

    debug!(
        "{} Longread BAM files, {} Shortread BAM files and {} assembly alignment BAMs {} Total BAMs",
        long_sample_count,
        short_sample_count,
        assembly_sample_count,
        (short_sample_count + long_sample_count + assembly_sample_count)
    );

    // Sliding window size for rolling SNV and SV counts
    // let window_size = m.value_of("window-size").unwrap().parse().unwrap();

    // The minimum contig size for binning
    let min_contig_size: u64 = m.value_of("min-contig-size").unwrap().parse().unwrap();

    // Set up multi progress bars
    // let multi = MultiProgress::new();
    let sty_eta = ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg} ETA: [{eta}]");

    let sty_aux = ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {spinner:.green} {msg} {pos:>4}/{len:4}");

    progress_bars[0].progress_bar.set_style(sty_eta.clone());

    let mut tree: Vec<&Elem> = Vec::with_capacity(progress_bars.len());

    for pb in progress_bars.iter() {
        tree.push(pb)
    }

    if (m.is_present("coverage-values")
        || (Path::new(&format!("{}/rosella_coverages.tsv", &output_prefix)).exists()
            && !m.is_present("force")))
        && (m.is_present("kmer-frequencies")
            || (Path::new(&format!("{}/rosella_kmer_table.tsv", &output_prefix)).exists()
                && !m.is_present("force")))
    {
        // Get the paths
        let mut coverage_path = None;
        if m.is_present("coverage-values") {
            coverage_path = Some(m.value_of("coverage-values").unwrap().to_string());
        } else if Path::new(&format!("{}/rosella_coverages.tsv", &output_prefix)).exists() {
            coverage_path = Some(format!("{}/rosella_coverages.tsv", &output_prefix));
        }

        // Get the paths
        let mut kmer_path = None;
        if m.is_present("kmer-frequencies") {
            kmer_path = Some(m.value_of("kmer-frequencies").unwrap().to_string());
        } else if Path::new(&format!("{}/rosella_kmer_table.tsv", &output_prefix)).exists() {
            kmer_path = Some(format!("{}/rosella_kmer_table.tsv", &output_prefix));
        }

        main_variant_matrix.read_inputs(coverage_path.as_deref(), kmer_path.as_deref(), None);

        // Completed contigs
        let pb = &progress_bars[0];
        // let pb = multi.insert(0, elem.progress_bar.clone());

        pb.progress_bar.enable_steady_tick(500);

        pb.progress_bar.finish_with_message(&format!(
            "Read results from previous run. If this is not desired please rerun with --force..."
        ));
        // multi.join().unwrap();
    } else if m.is_present("coverage-values") && m.is_present("kmer-frequencies") {
        // Read from provided inputs
        main_variant_matrix.read_inputs(
            Some(m.value_of("coverage-values").unwrap()),
            Some(m.value_of("kmer-frequencies").unwrap()),
            // Some(m.value_of("variant-rates").unwrap()),
            None,
        );
        // Completed contigs
        let pb = &progress_bars[0];
        // let pb = multi.insert(0, elem.progress_bar.clone());

        pb.progress_bar.enable_steady_tick(500);

        pb.progress_bar.finish_with_message(&format!("Using provided input files..."));
        // multi.join().unwrap();
    } else if short_sample_count + long_sample_count > 0 {
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
        let contig_header = match finish_bams(bam_readers, n_threads) {
            Some(header) => header,
            None => panic!("No header retrieved from BAM files"),
        };
        let n_contigs = contig_header.target_count();

        let contig_lens: HashMap<u32, u64> =
            (0..n_contigs)
                .into_iter()
                .fold(HashMap::new(), |mut map: HashMap<u32, u64>, tid| {
                    map.insert(tid, contig_header.target_len(tid).unwrap());
                    map
                });

        let contig_names = contig_header.target_names();
        {
            // Completed contigs
            let pb = &progress_bars[0];
            // let pb = multi.insert(0, elem.progress_bar.clone());

            pb.progress_bar.enable_steady_tick(500);

            pb.progress_bar.set_message(&format!("{}...", &pb.key,));
        }

        let mut indexed_reference = generate_faidx(reference);

        for tid in (0..n_contigs).into_iter() {
            let tmp_bam_file_cache = match tmp_bam_file_cache.as_ref() {
                Some(cache) => Some(cache.path().to_str().unwrap().to_string()),
                None => None,
            };

            let target_name = str::from_utf8(contig_names[tid as usize]).unwrap();

            if contig_lens.get(&tid).unwrap() < &1000 {
                {
                    let pb = &progress_bars;

                    pb[0].progress_bar.inc(1);
                    pb[0]
                        .progress_bar
                        .set_message(&format!("{} analyzed...", target_name));
                    let pos = pb[0].progress_bar.position();
                    let len = pb[0].progress_bar.length();
                    if pos >= len {
                        pb[0]
                            .progress_bar
                            .finish_with_message(&format!("All genomes analyzed {}", "✔",));
                    }
                }
                continue;
            }

            // let mut coverage_estimators = coverage_estimators.clone();

            let mut per_reference_samples =
                short_sample_count + long_sample_count + assembly_sample_count;
            let mut per_reference_short_samples = short_sample_count;
            let mut variant_matrix = VariantMatrix::new_matrix(
                short_sample_count + long_sample_count + assembly_sample_count,
                None,
            );

            // // Read BAMs back in as indexed
            let mut indexed_bam_readers = recover_bams(
                m,
                short_sample_count,
                long_sample_count,
                assembly_sample_count,
                &tmp_bam_file_cache,
            );

            // let mut variant_matrix = Mutex::new(variant_matrix);
            if !m.is_present("coverage-values") {
                // We just grab the coverage
                indexed_bam_readers.into_iter().enumerate().for_each(
                    |(sample_idx, bam_generator)| {
                        let mut bam_generator = generate_indexed_named_bam_readers_from_bam_files(
                            vec![&bam_generator],
                            n_threads as u32,
                        )
                        .into_iter()
                        .next()
                        .unwrap();
                        if sample_idx < short_sample_count {
                            contig_coverage(
                                bam_generator,
                                sample_idx,
                                per_reference_samples,
                                coverage_estimators,
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
                                tid,
                            )
                        } else if sample_idx >= short_sample_count
                            && sample_idx < (short_sample_count + long_sample_count)
                        {
                            contig_coverage(
                                bam_generator,
                                sample_idx,
                                per_reference_samples,
                                coverage_estimators,
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
                                tid,
                            )
                        } else if sample_idx >= (short_sample_count + long_sample_count) {
                            // Skip assembly bams here
                        }
                    },
                );
            }

            // Collects info about variants across samples to check whether they are genuine or not
            // using FDR

            // variant_matrix.remove_false_discoveries(tid as i32, alpha, reference);

            // K-mer size for kmer frequency table
            let kmer_size: usize = m.value_of("kmer-size").unwrap().parse().unwrap();
            {
                main_variant_matrix.merge_matrices(
                    tid,
                    variant_matrix,
                    (short_sample_count + long_sample_count),
                );

                if !m.is_present("kmer-frequencies") {
                    main_variant_matrix.calc_kmer_frequencies(
                        tid,
                        kmer_size,
                        &mut indexed_reference,
                        n_contigs as usize,
                    );
                }
            }
            {
                let pb = &progress_bars;

                pb[0].progress_bar.inc(1);
                pb[0]
                    .progress_bar
                    .set_message(&format!("{} analyzed...", target_name));
                let pos = pb[0].progress_bar.position();
                let len = pb[0].progress_bar.length();
                if pos >= len {
                    pb[0]
                        .progress_bar
                        .finish_with_message(&format!("All genomes analyzed {}", "✔",));
                }
            }
        }
        // multi.join().unwrap();
    } else if m.is_present("coverage-values") {
        // Replace coverage values with CoverM values if available
        main_variant_matrix.read_inputs(Some(m.value_of("coverage-values").unwrap()), None, None);
        let n_contigs = main_variant_matrix.get_n_contigs();
        let contig_lens = main_variant_matrix.get_contig_lengths();

        {
            // Completed contigs
            let pb = &progress_bars[0];
            // let pb = multi.insert(0, elem.progress_bar.clone());

            pb.progress_bar.enable_steady_tick(500);

            pb.progress_bar.set_message(&format!("{}...", &pb.key,));
        }

        let mut indexed_reference = generate_faidx(reference);

        for tid in (0..n_contigs).into_iter() {
            if contig_lens.get(&(tid as i32)).unwrap() < &1000 {
                {
                    let pb = &tree;

                    pb[0].progress_bar.inc(1);
                    pb[0].progress_bar.set_message(&format!("analyzed..."));
                    let pos = pb[0].progress_bar.position();
                    let len = pb[0].progress_bar.length();
                    if pos >= len {
                        pb[0]
                            .progress_bar
                            .finish_with_message(&format!("All contigs analyzed {}", "✔",));
                    }
                }
                continue;
            }
            {
                if !m.is_present("kmer-frequencies") {
                    // K-mer size for kmer frequency table
                    let kmer_size: usize = m.value_of("kmer-size").unwrap().parse().unwrap();
                    main_variant_matrix.calc_kmer_frequencies(
                        tid,
                        kmer_size,
                        &mut indexed_reference,
                        n_contigs as usize,
                    );
                }
            }
            {
                let pb = &progress_bars;

                pb[0].progress_bar.inc(1);
                pb[0].progress_bar.set_message(&format!("analyzed..."));
                let pos = pb[0].progress_bar.position();
                let len = pb[0].progress_bar.length();
                if pos >= len {
                    pb[0]
                        .progress_bar
                        .finish_with_message(&format!("All contigs analyzed {}", "✔",));
                }
            }
        }
        // multi.join().unwrap();
    } else {
        warn!(
            "ERROR: User has not supplied reads, BAM files, coverage results \
        or directory containing results of previous run. Cannot proceed. Exiting..."
        );
        process::exit(1)
    }

    let pb = Elem {
        key: "Binning",
        index: 0,
        progress_bar: ProgressBar::new(6),
    };

    pb.progress_bar.set_style(sty_aux.clone());
    pb.progress_bar.enable_steady_tick(500);
    std::fs::create_dir_all(&output_prefix).unwrap();

    for i in (0..3).into_iter() {
        if i == 0 {
            if !m.is_present("kmer-frequences") {
                if !Path::new(&format!("{}/rosella_kmer_table.tsv", &output_prefix)).exists()
                    || m.is_present("force")
                {
                    main_variant_matrix.write_kmer_table(&output_prefix, min_contig_size);
                }
            }
            pb.progress_bar.inc(1);
            pb.progress_bar
                .set_message(&format!("K-mer frequencies written {}", "✔",));
        } else if i == 1 {
            // if !m.is_present("variant-rates") {
            //     if !Path::new(&format!("{}/rosella_variant_rates.tsv", &output_prefix))
            //         .exists()
            //         || m.is_present("force")
            //     {
            //         main_variant_matrix.write_variant_rates(&output_prefix);
            //     }
            // }
            {
                pb.progress_bar.inc(1);
                pb.progress_bar
                    .set_message(&format!("Variant rates written {}", "✔",));
            }
        } else if i == 2 {
            if !m.is_present("coverage-values") {
                if !Path::new(&format!("{}/rosella_coverages.tsv", &output_prefix)).exists()
                    || m.is_present("force")
                {
                    main_variant_matrix.write_coverage(&output_prefix);
                }
            }
            {
                pb.progress_bar.inc(1);
                pb.progress_bar
                    .set_message(&format!("Contig coverages written {}", "✔",));
            }
        }
    }

    pb.progress_bar
        .set_message(&format!("Calculating UMAP embeddings and clustering...",));
    main_variant_matrix.bin_contigs(output_prefix, m);
    pb.progress_bar.inc(1);

    pb.progress_bar.set_message(&format!("Finalizing bins...",));
    main_variant_matrix.finalize_bins(output_prefix, m);
    pb.progress_bar.inc(1);

    pb.progress_bar.set_message(&format!("Writing bins...",));
    main_variant_matrix.write_bins(output_prefix, &reference);
    pb.progress_bar.inc(1);

    pb.progress_bar
        .finish_with_message(&format!("All steps completed {}", "✔",));

    info!("Analysis finished!");
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
