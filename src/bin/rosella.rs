extern crate openssl;
extern crate openssl_sys;

extern crate rosella;

use rosella::cli::*;
use rosella::estimation::contig;
use rosella::external_command_checker;
use rosella::utils::*;
use rosella::*;

extern crate rust_htslib;

extern crate bio;
use bio::alignment::sparse::*;

extern crate coverm;
use coverm::bam_generator::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::FlagFilter;
use coverm::*;

use std::collections::BTreeMap;
use std::env;
use std::path::Path;
use std::process;
use std::str;

extern crate tempdir;
extern crate tempfile;
use tempfile::NamedTempFile;
extern crate clap;
use clap::*;

extern crate itertools;
use itertools::Itertools;

#[macro_use]
extern crate log;
use log::LevelFilter;
extern crate env_logger;
use env_logger::Builder;

fn main() {
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    set_log_level(&matches, false);

    match matches.subcommand_name() {
        Some("recover") => {
            let m = matches.subcommand_matches("recover").unwrap();
            let mode = "bin";
            if m.is_present("full-help") {
                println!("{}", binning_full_help());
                process::exit(1);
            }
            prepare_pileup(m, mode);
        },
        Some("refine") => {
            let m = matches.subcommand_matches("refine").unwrap();
            let mode = "refine";
            if m.is_present("full-help") {
                println!("{}", refining_full_help());
                process::exit(1);
            }
            prepare_pileup(m, mode);
        }
        Some("kmer") => {
            let m = matches.subcommand_matches("kmer").unwrap();
            if m.is_present("full-help") {
                //                println!("{}", contig_full_help());
                process::exit(1);
            }
            set_log_level(m, true);
            let reference_path = Path::new(m.value_of("reference").unwrap());
            let fasta_reader = match bio::io::fasta::Reader::from_file(&reference_path) {
                Ok(reader) => reader,
                Err(e) => {
                    eprintln!("Missing or corrupt fasta file {}", e);
                    process::exit(1);
                }
            };
            let contigs = fasta_reader.records().collect_vec();
            // Initialize bound contig variable
            let mut tet_freq = BTreeMap::new();
            let contig_count = contigs.len();
            let mut contig_idx = 0 as usize;
            let mut contig_names = vec![String::new(); contig_count];
            for contig in contigs {
                let contig = contig.unwrap();
                contig_names[contig_idx] = contig.id().to_string();
                debug!("Parsing contig: {}", contig.id());
                let kmers = hash_kmers(contig.seq(), 4);
                // Get kmer counts in a contig
                for (kmer, pos) in kmers {
                    let k = tet_freq
                        .entry(kmer.to_vec())
                        .or_insert(vec![0; contig_count]);
                    k[contig_idx] = pos.len();
                }
                contig_idx += 1;
            }
            for (idx, contig) in contig_names.iter().enumerate() {
                print!("{}\t", contig);
                for (_kmer, counts) in &tet_freq {
                    print!("{}\t", counts[idx])
                }
                print!("\n");
            }
        }
        _ => {
            app.print_help().unwrap();
            println!();
        }
    }
}

fn prepare_pileup(m: &clap::ArgMatches, mode: &str) {
    // This function is amazingly painful. It handles every combination of longread and short read
    // mapping or bam file reading. Could not make it smaller using dynamic or static dispatch
    set_log_level(m, true);
    let estimators = EstimatorsAndTaker::generate_from_clap(m);
    let filter_params = FilterParameters::generate_from_clap(m);
    let threads = m.value_of("threads").unwrap().parse().unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    prepare_binning(
        m,
        mode,
        estimators,
        filter_params,
        threads
    );

}

fn prepare_binning(
    m: &clap::ArgMatches,
    mode: &str,
    mut estimators: EstimatorsAndTaker,
    filter_params: FilterParameters,
    _threads: usize,
) {
    // Temp directory that will house all cached bams for variant calling
    let tmp_dir = match m.is_present("bam-file-cache-directory") {
        false => {
            let tmp_direct = tempdir::TempDir::new("rosella_fifo")
                .expect("Unable to create temporary directory");
            debug!("Temp directory {}", tmp_direct.as_ref().to_str().unwrap());
            std::fs::create_dir(format!("{}/long", &tmp_direct.as_ref().to_str().unwrap()))
                .unwrap();
            std::fs::create_dir(format!("{}/short", &tmp_direct.as_ref().to_str().unwrap()))
                .unwrap();
            std::fs::create_dir(format!(
                "{}/assembly",
                &tmp_direct.as_ref().to_str().unwrap()
            ))
                .unwrap();

            Some(tmp_direct)
        }
        true => None,
    };

    // let (concatenated_genomes, genomes_and_contigs_option) = setup_genome_fasta_files(&m);
    // debug!("Found genomes_and_contigs {:?}", genomes_and_contigs_option);
    let references = match mode {
        "bin" => vec![m.value_of("reference").unwrap()],
        "refine" => vec![m.value_of("assembly").unwrap()],
        _ => vec![m.value_of("reference").unwrap()]
    };
    if m.is_present("bam-files") {
        let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();

        // Associate genomes and contig names, if required
        if filter_params.doing_filtering() {
            let bam_readers = bam_generator::generate_filtered_bam_readers_from_bam_files(
                bam_files,
                filter_params.flag_filters.clone(),
                filter_params.min_aligned_length_single,
                filter_params.min_percent_identity_single,
                filter_params.min_aligned_percent_single,
                filter_params.min_aligned_length_pair,
                filter_params.min_percent_identity_pair,
                filter_params.min_aligned_percent_pair,
            );

            if m.is_present("longread-bam-files") {
                let bam_files = m.values_of("longread-bam-files").unwrap().collect();
                let long_readers =
                    bam_generator::generate_named_bam_readers_from_bam_files(bam_files);

                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    Some(bam_readers),
                    filter_params.flag_filters,
                    Some(long_readers),
                    None::<Vec<PlaceholderBamFileReader>>,
                    None,
                    tmp_dir,
                    None,
                )
            } else if m.is_present("longreads") {
                // Perform mapping
                let long_generators =
                    long_generator_setup(&m, &None, &Some(references.clone()), &tmp_dir);

                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    Some(bam_readers),
                    filter_params.flag_filters,
                    Some(long_generators),
                    None::<Vec<PlaceholderBamFileReader>>,
                    None,
                    tmp_dir,
                    None,
                )
            } else {
                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    Some(bam_readers),
                    filter_params.flag_filters,
                    None::<Vec<PlaceholderBamFileReader>>,
                    None::<Vec<PlaceholderBamFileReader>>,
                    None,
                    tmp_dir,
                    None,
                )
            }
        } else {
            let bam_readers = bam_generator::generate_named_bam_readers_from_bam_files(bam_files);

            if m.is_present("longread-bam-files") {
                let bam_files = m.values_of("longread-bam-files").unwrap().collect();
                let long_readers =
                    bam_generator::generate_named_bam_readers_from_bam_files(bam_files);
                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    Some(bam_readers),
                    filter_params.flag_filters,
                    Some(long_readers),
                    None::<Vec<PlaceholderBamFileReader>>,
                    None,
                    tmp_dir,
                    None,
                )
            } else if m.is_present("longreads") {
                // Perform mapping
                let long_generators =
                    long_generator_setup(&m, &None, &Some(references.clone()), &tmp_dir);

                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    Some(bam_readers),
                    filter_params.flag_filters,
                    Some(long_generators),
                    None::<Vec<PlaceholderBamFileReader>>,
                    None,
                    tmp_dir,
                    None,
                )
            } else {
                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    Some(bam_readers),
                    filter_params.flag_filters,
                    None::<Vec<PlaceholderBamFileReader>>,
                    None::<Vec<PlaceholderBamFileReader>>,
                    None,
                    tmp_dir,
                    None,
                )
            }
        }
    } else if m.is_present("read1")
        || m.is_present("coupled")
        || m.is_present("interleaved")
        || m.is_present("single")
    {
        let mapping_program = parse_mapping_program(m.value_of("mapper"));
        external_command_checker::check_for_samtools();

        if filter_params.doing_filtering() {
            debug!("Filtering..");
            let readtype = ReadType::Short;
            let generator_sets = get_streamed_filtered_bam_readers(
                m,
                mapping_program,
                &None,
                &filter_params,
                &readtype,
                &Some(references.clone()),
                &tmp_dir,
            );
            let mut all_generators = vec![];
            let mut indices = vec![]; // Prevent indices from being dropped
            for set in generator_sets {
                indices.push(set.index);
                for g in set.generators {
                    all_generators.push(g)
                }
            }
            debug!("Finished collecting generators.");
            if m.is_present("longread-bam-files") {
                let bam_files = m.values_of("longread-bam-files").unwrap().collect();
                let long_readers =
                    bam_generator::generate_named_bam_readers_from_bam_files(bam_files);

                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    Some(all_generators),
                    filter_params.flag_filters,
                    Some(long_readers),
                    None::<Vec<PlaceholderBamFileReader>>,
                    None,
                    tmp_dir,
                    None,
                )
            } else if m.is_present("longreads") {
                // Perform mapping
                let long_generators =
                    long_generator_setup(&m, &None, &Some(references.clone()), &tmp_dir);
                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    Some(all_generators),
                    filter_params.flag_filters,
                    Some(long_generators),
                    None::<Vec<PlaceholderBamFileReader>>,
                    None,
                    tmp_dir,
                    None,
                )
            } else {
                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    Some(all_generators),
                    filter_params.flag_filters,
                    None::<Vec<PlaceholderBamFileReader>>,
                    None::<Vec<PlaceholderBamFileReader>>,
                    None,
                    tmp_dir,
                    None,
                )
            }
        } else {
            debug!("Not filtering..");
            let readtype = ReadType::Short;
            let generator_sets = get_streamed_bam_readers(
                m,
                mapping_program,
                &None,
                &readtype,
                &Some(references.clone()),
                &tmp_dir,
            );
            let mut all_generators = vec![];
            let mut indices = vec![]; // Prevent indices from being dropped
            for set in generator_sets {
                indices.push(set.index);
                for g in set.generators {
                    all_generators.push(g)
                }
            }
            if m.is_present("longread-bam-files") {
                let bam_files = m.values_of("longread-bam-files").unwrap().collect();
                let long_readers =
                    bam_generator::generate_named_bam_readers_from_bam_files(bam_files);

                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    Some(all_generators),
                    filter_params.flag_filters,
                    Some(long_readers),
                    None::<Vec<PlaceholderBamFileReader>>,
                    None,
                    tmp_dir,
                    None,
                )
            } else if m.is_present("longreads") {
                // Perform mapping
                let long_generators =
                    long_generator_setup(&m, &None, &Some(references.clone()), &tmp_dir);

                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    Some(all_generators),
                    filter_params.flag_filters,
                    Some(long_generators),
                    None::<Vec<PlaceholderBamFileReader>>,
                    None,
                    tmp_dir,
                    None,
                )
            } else {
                run_pileup(
                    m,
                    mode,
                    &mut estimators,
                    Some(all_generators),
                    filter_params.flag_filters,
                    None::<Vec<PlaceholderBamFileReader>>,
                    None::<Vec<PlaceholderBamFileReader>>,
                    None,
                    tmp_dir,
                    None,
                )
            }
        }
    } else {
        run_pileup(
            m,
            mode,
            &mut estimators,
            None::<Vec<PlaceholderBamFileReader>>,
            filter_params.flag_filters,
            None::<Vec<PlaceholderBamFileReader>>,
            None::<Vec<PlaceholderBamFileReader>>,
            None,
            tmp_dir,
            None,
        )
    }
}

struct EstimatorsAndTaker {
    estimators: Vec<CoverageEstimator>,
}

impl EstimatorsAndTaker {
    pub fn generate_from_clap(m: &clap::ArgMatches) -> EstimatorsAndTaker {
        let mut estimators = vec![];
        let min_fraction_covered = parse_percentage(&m, "min-covered-fraction");
        let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();

        let methods: Vec<&str> = m.values_of("method").unwrap().collect();

        if doing_metabat(&m) {
            estimators.push(CoverageEstimator::new_estimator_length());
            estimators.push(CoverageEstimator::new_estimator_mean(
                min_fraction_covered,
                contig_end_exclusion,
                false,
            ));
            estimators.push(CoverageEstimator::new_estimator_variance(
                min_fraction_covered,
                contig_end_exclusion,
            ));

            debug!("Cached regular coverage taker for metabat mode being used");
        } else {
            for (_i, method) in methods.iter().enumerate() {
                match method {
                    &"mean" => {
                        estimators.push(CoverageEstimator::new_estimator_length());

                        estimators.push(CoverageEstimator::new_estimator_mean(
                            min_fraction_covered,
                            contig_end_exclusion,
                            false,
                        )); // TODO: Parameterise exclude_mismatches

                        estimators.push(CoverageEstimator::new_estimator_variance(
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));
                    }
                    &"trimmed_mean" => {
                        let min = value_t!(m.value_of("trim-min"), f32).unwrap();
                        let max = value_t!(m.value_of("trim-max"), f32).unwrap();
                        if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                            error!(
                                "error: Trim bounds must be between 0 and 1, and \
                                 min must be less than max, found {} and {}",
                                min, max
                            );
                            process::exit(1);
                        }
                        estimators.push(CoverageEstimator::new_estimator_length());

                        estimators.push(CoverageEstimator::new_estimator_trimmed_mean(
                            min,
                            max,
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));

                        estimators.push(CoverageEstimator::new_estimator_variance(
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));
                    }
                    _ => unreachable!(),
                };
            }
        }

        // Check that min-covered-fraction is being used as expected
        if min_fraction_covered != 0.0 {
            let die = |estimator_name| {
                error!(
                    "The '{}' coverage estimator cannot be used when \
                     --min-covered-fraction is > 0 as it does not calculate \
                     the covered fraction. You may wish to set the \
                     --min-covered-fraction to 0 and/or run this estimator \
                     separately.",
                    estimator_name
                );
                process::exit(1)
            };
            for e in &estimators {
                match e {
                    CoverageEstimator::ReadCountCalculator { .. } => die("counts"),
                    CoverageEstimator::ReferenceLengthCalculator { .. } => die("length"),
                    CoverageEstimator::ReadsPerBaseCalculator { .. } => die("reads_per_base"),
                    _ => {}
                }
            }
        }

        return EstimatorsAndTaker {
            estimators: estimators,
        };
    }
}

fn run_pileup<
    'a,
    R: NamedBamReader,
    S: NamedBamReaderGenerator<R>,
    T: NamedBamReader,
    U: NamedBamReaderGenerator<T>,
    V: NamedBamReader,
    W: NamedBamReaderGenerator<V>,
>(
    m: &clap::ArgMatches,
    mode: &str,
    estimators: &mut EstimatorsAndTaker,
    bam_readers: Option<Vec<S>>,
    flag_filters: FlagFilter,
    long_readers: Option<Vec<U>>,
    assembly_readers: Option<Vec<W>>,
    _genomes_and_contigs_option: Option<GenomesAndContigs>,
    tmp_bam_file_cache: Option<tempdir::TempDir>,
    _concatenated_genomes: Option<NamedTempFile>,
) {
    match mode {
        "bin" | "refine" => {
            let mapq_threshold = m.value_of("mapq-threshold").unwrap().parse().unwrap();

            let output_prefix = match m.is_present("output-directory") {
                true => {
                    match std::fs::create_dir_all(
                        m.value_of("output-directory").unwrap().to_string(),
                    ) {
                        Ok(_) => {}
                        Err(err) => panic!("Unable to create output directory {:?}", err),
                    };
                    m.value_of("output-directory").unwrap()
                }
                false => "./",
            };

            let threads = m.value_of("threads").unwrap().parse().unwrap();

            let method = m.value_of("method").unwrap();

            let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();
            let min = value_t!(m.value_of("trim-min"), f32).unwrap();
            let max = value_t!(m.value_of("trim-max"), f32).unwrap();
            if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                panic!(
                    "error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}",
                    min, max
                );
            }

            contig::prepare_flight(
                m,
                bam_readers,
                long_readers,
                assembly_readers,
                mode,
                &mut estimators.estimators,
                flag_filters,
                mapq_threshold,
                min,
                max,
                contig_end_exclusion,
                output_prefix,
                threads,
                method,
                m.is_present("longread-bam-files"),
                tmp_bam_file_cache,
            );
        }
        _ => panic!("Unknown rosella mode"),
    }
}

fn set_log_level(matches: &clap::ArgMatches, is_last: bool) {
    let mut log_level = LevelFilter::Info;
    let mut specified = false;
    if matches.is_present("verbose") {
        specified = true;
        log_level = LevelFilter::Debug;
    }
    if matches.is_present("quiet") {
        specified = true;
        log_level = LevelFilter::Error;
    }
    if specified || is_last {
        let mut builder = Builder::new();
        builder.filter_level(log_level);
        if env::var("RUST_LOG").is_ok() {
            builder.parse_filters(&env::var("RUST_LOG").unwrap());
        }
        if builder.try_init().is_err() {
            panic!("Failed to set log level - has it been specified multiple times?")
        }
    }
    if is_last {
        info!("rosella version {}", crate_version!());
    }
}
