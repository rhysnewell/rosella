use bird_tool_utils::clap_utils::{default_roff, monospace_roff, table_roff};
use bird_tool_utils_man::prelude::{Author, Flag, Manual, Opt, Section};
use clap::crate_version;

fn add_mapping_options(manual: Manual) -> Manual {
    manual.custom(
        Section::new("Mapping algorithm options")
            .option(Opt::new("NAME").long("--mapper").help(&format!(
                "Underlying mapping software used for short reads {}. One of: {}",
                default_roff("minimap2-sr"),
                table_roff(&[
                    &["name", "description"],
                    &[
                        &monospace_roff("minimap2-sr"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x sr"))
                    ],
                    &[
                        &monospace_roff("bwa-mem"),
                        &format!("bwa mem using default parameters")
                    ],
                    &[
                        &monospace_roff("bwa-mem2"),
                        &format!("bwa-mem2 using default parameters")
                    ],
                    &[
                        &monospace_roff("minimap2-ont"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-ont"))
                    ],
                    &[
                        &monospace_roff("minimap2-pb"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-pb"))
                    ],
                    &[
                        &monospace_roff("minimap2-hifi"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-hifi"))
                    ],
                    &[
                        &monospace_roff("minimap2-no-preset"),
                        &format!("minimap2 with no '{}' option", &monospace_roff("-x"))
                    ],
                ])
            )))
            .option(Opt::new("NAME").long("--longread-mapper").help(&format!(
                "Underlying mapping software used for long reads {}. One of: {}",
                default_roff("minimap2-ont"),
                table_roff(&[
                    &["name", "description"],
                    &[
                        &monospace_roff("minimap2-ont"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-ont"))
                    ],
                    &[
                        &monospace_roff("minimap2-pb"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-pb"))
                    ],
                    &[
                        &monospace_roff("minimap2-hifi"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-hifi"))
                    ],
                    &[
                        &monospace_roff("minimap2-no-preset"),
                        &format!("minimap2 with no '{}' option", &monospace_roff("-x"))
                    ],
                ])
            )))
            .option(Opt::new("PARAMS").long("--minimap2-params").help(&format!(
                "Extra parameters to provide to minimap2, \
        both indexing command (if used) and for \
        mapping. Note that usage of this parameter \
        has security implications if untrusted input \
        is specified. '{}' is always specified to minimap2. \
        [default: none] \n",
                &monospace_roff("-a")
            )))
            .flag(Flag::new().long("--minimap2-reference-is-index").help(
                "Treat reference as a minimap2 database, not as a FASTA file. [default: not set]",
            ))
            .option(Opt::new("PARAMS").long("--bwa-params").help(
                "Extra parameters to provide to BWA or BWA-MEM2. Note \
        that usage of this parameter has security \
        implications if untrusted input is specified. \
        [default: none] \n",
            )),
    )
}

fn add_thresholding_options(manual: Manual) -> Manual {
    manual.custom(
        Section::new("Alignment thresholding")
            .option(
                Opt::new("INT")
                    .long("--min-read-aligned-length")
                    .help(&format!(
                        "Exclude reads with smaller numbers of \
        aligned bases. {} \n",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-percent-identity")
                    .help(&format!(
                        "Exclude reads by overall percent \
        identity e.g. 95 for 95%. {} \n",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-aligned-percent")
                    .help(&format!(
                        "Exclude reads by percent aligned \
        bases e.g. 95 means 95% of the read's \
        bases must be aligned. {} \n",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("INT")
                    .long("--min-read-aligned-length-pair")
                    .help(&format!(
                        "Exclude pairs with smaller numbers of \
        aligned bases. \
        Implies --proper-pairs-only. {} \n",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-percent-identity-pair")
                    .help(&format!(
                        "Exclude pairs by overall percent \
                identity e.g. 95 for 95%. \
                Implies --proper-pairs-only. {} \n",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-aligned-percent-pair")
                    .help(&format!(
                        "Exclude reads by percent aligned \
                bases e.g. 95 means 95% of the read's \
                bases must be aligned. \
                Implies --proper-pairs-only. {} \n",
                        default_roff("0")
                    )),
            )
            .flag(
                Flag::new()
                    .long("--proper-pairs-only")
                    .help("Require reads to be mapped as proper pairs. [default: not set] \n"),
            )
            .flag(
                Flag::new()
                    .long("--exclude-supplementary")
                    .help("Exclude supplementary alignments. [default: not set] \n"),
            )
            .flag(
                Flag::new()
                    .long("--include-secondary")
                    .help("Include secondary alignments. [default: not set] \n"),
            )
            .option(Opt::new("INT").long("--contig-end-exclusion").help(
                "Exclude bases at the ends of reference \n
                         sequences from calculation [default: 75]",
            ))
            .option(Opt::new("FLOAT").long("--trim-min").help(
                "Remove this smallest fraction of positions \n
                         when calculating trimmed_mean [default: 5]",
            ))
            .option(Opt::new("FLOAT").long("--trim-max").help(
                "Maximum fraction for trimmed_mean \n
                         calculations [default: 95]",
            ))
    )
}

fn read_mapping_params_section() -> Section {
    Section::new("Read mapping parameters")
        .option(
            Opt::new("PATH ..")
                .short("-1")
                .help("Forward FASTA/Q file(s) for mapping. These may be gzipped or not. \n"),
        )
        .option(
            Opt::new("PATH ..")
                .short("-2")
                .help("Reverse FASTA/Q file(s) for mapping. These may be gzipped or not. \n"),
        )
        .option(Opt::new("PATH ..").short("-c").long("--coupled").help(
            "One or more pairs of forward and reverse \
        possibly gzipped FASTA/Q files for mapping in order \
        <sample1_R1.fq.gz> <sample1_R2.fq.gz> \
        <sample2_R1.fq.gz> <sample2_R2.fq.gz> .. \n",
        ))
        .option(
            Opt::new("PATH ..")
                .long("--interleaved")
                .help("Interleaved FASTA/Q files(s) for mapping. These may be gzipped or not. \n"),
        )
        .option(
            Opt::new("PATH ..")
                .long("--single")
                .help("Unpaired FASTA/Q files(s) for mapping. These may be gzipped or not. \n"),
        )
        .option(
            Opt::new("PATH ..")
                .long("--longreads")
                .help("Longread FASTA/Q files(s) for mapping. These may be gzipped or not. \n"),
        )
        .option(
            Opt::new("PATH")
                .short("-b")
                .long("--bam-files")
                .help(&format!(
                    "Path to BAM file(s). These must be \
                reference sorted (e.g. with samtools sort) \
                unless {} is specified, in which \
                case they must be read name sorted (e.g. \
                with {}). When specified, no read mapping algorithm is undertaken. \n",
                    monospace_roff("--sharded"),
                    monospace_roff("samtools sort -n"),
                )),
        )
        .option(
            Opt::new("PATH")
                .short("-l")
                .long("--longread-bam-files")
                .help(&format!(
                    "Path to longread BAM file(s). These must be \
                reference sorted (e.g. with samtools sort) \
                unless {} is specified, in which \
                case they must be read name sorted (e.g. \
                with {}). When specified, no read mapping algorithm is undertaken. \n",
                    monospace_roff("--sharded"),
                    monospace_roff("samtools sort -n"),
                )),
        )
}

fn binning_params_section() -> Section {
    Section::new("Binning parameters")
        .option(
            Opt::new("PATH")
                .short("-o")
                .long("--output-directory")
                .help(&format!(
                    "Output directory for binning results. \
                [default: {}] \n",
                    monospace_roff("rosella_output")
                )),
        )
        .option(
            Opt::new("PATH")
                .short("-C")
                .long("--coverage-file")
                .help(&format!(
                    "The output from the results of {} \
                in MetaBAT mode on the provided assembly and short read samples. \
                If not provided, rosella will calculate coverage values. \n",
                    monospace_roff("CoverM contig")
                )),
        )
        .option(Opt::new("PATH").long("--kmer-frequency-file")
            .short("-K")
            .help(&format!(
                "The kmer frequency table created by {}. \n",
                monospace_roff("rosella")
            ))
        )
        .option(Opt::new("INT").long("--min-contig-size").help(&format!(
            "Minimum contig size in base pairs to be considered for binning. \
                Contigs between 1000 bp and this value will be recovered in \
                the contig rescue stage if multiple samples are available. \
                [default: {}] \n",
            default_roff("1500")
        )))
        .option(Opt::new("INT").long("--min-bin-size").help(&format!(
            "Minimum bin size in base pairs for MAG to be reported. If a bin \
                is smaller than this, then it will be split apart and the contigs \
                will be added to other bins. [default: {}] \n",
            default_roff("100000")
        )))
        .flag(
            Flag::new()
                .long("--refine")
                .help("Re-cluster each recovered bin on its own and keep any split that \
                    leaves the pieces tighter than the whole. Raises purity at the cost \
                    of completeness, so it is off by default. \n"),
        )
        .option(Opt::new("INT").long("--max-bin-size").help(&format!(
            "Bin size in base pairs above which a bin is always re-clustered \
                during refinement, however good its internal statistics look. \
                [default: {}] \n",
            default_roff("15000000")
        )))
        .option(
            Opt::new("INT")
                .long("--kmer-size")
                .short("-k")
                .help(&format!(
                    "Kmer size to use for kmer frequency table. [default: {}] \n",
                    default_roff("4")
                )),
        )
        .option(
            Opt::new("--n-neighbours")
                .long("--n-neighbours")
                .short("-n")
                .help(&format!(
                    "Number of neighbors used in the UMAP algorithm. \
                [default: {}] \n",
                    default_roff("100")
                )),
        )
        .option(
            Opt::new("INT")
                .long("--max-retries")
                .help(&format!(
                    "Maximum number of times to retry refining a genome \
                    if it fails. [default: {}] \n",
                    default_roff("5")
                )),
        )
}

fn refining_options() -> Section {
    Section::new("Genome refining options")
        .option(
            Opt::new("PATH")
                .short("-d")
                .long("--genome-fasta-directory")
                .help(&format!(
                    "Directory containing FASTA files of contigs \
                    [required unless {} is specified] \n",
                    monospace_roff("-f/--genome-fasta-files")
                )),
        )
        .option(
            Opt::new("PATH")
                .short("-f")
                .long("--genome-fasta-files")
                .help(&format!(
                    "FASTA files of contigs e.g. concatenated \
                    genomes or metagenome assembly
                    [required unless {} is specified] \n",
                    monospace_roff("-i/--genome-fasta-directory")
                )),
        )
        .option(
            Opt::new("STR")
                .short("-x")
                .long("--genome-fasta-extension")
                .help(&format!(
                    "FASTA file extension in --genome-fasta-directory \
                        [default \"fna\"] \n"
                )),
        )
        .option(
            Opt::new("PATH")
                .long("--checkm-results")
                .help(&format!(
                    "CheckM 1 or 2 results that contain information on the \
                    completeness and contamination of the input genomes. \n"
                )),
        )
        .option(
            Opt::new("FLOAT")
                .long("--max-contamination")
                .help(&format!(
                    "Contamination threshold used when CheckM results are \
                    provided. Genomes above it are deconstructed and rebuilt, \
                    genomes at or below it are left unchanged. [default: {}] \n",
                    default_roff("15.0")
                )),
        )
        .option(
            Opt::new("STR")
                .long("--bin-tag")
                .help(&format!(
                    "Tag to use for the refined bins. [default: {}] \n",
                    default_roff("refined_1")
                )),
        )
}

fn reference_options_simple() -> Section {
    Section::new("Input assembly option")
        .option(
            Opt::new("PATH")
                .short("-r,-f")
                .long("--assembly,--reference")
                .help(&format!(
                    "FASTA files of contigs e.g. concatenated \
                    genomes or metagenome assembly
                    [required unless {} is specified] \n",
                    monospace_roff("-d/--genome-fasta-directory")
                )),
        )
}

fn threads_options() -> Section {
    Section::new("Threading options").option(
        Opt::new("INT")
            .long("--threads")
            .short("-t")
            .help("Maximum number of threads used. [default: 10] \n"),
    )
}

fn add_verbosity_flags(manual: Manual) -> Manual {
    manual
        .flag(
            Flag::new()
                .short("-v")
                .long("--verbose")
                .help("Print extra debugging information. [default: not set] \n"),
        )
        .flag(Flag::new().short("-q").long("--quiet").help(
            "Unless there is an error, do not print \
    log messages. [default: not set]",
        ))
}

pub fn recover_full_help() -> Manual {
    let mut manual = Manual::new("rosella recover")
        .about(
            &format!(
                "Recover MAGs from contigs using UMAP and HDBSCAN clustering. (version {})",
                crate_version!()
            )
        )
        .author(Author::new(crate::AUTHOR).email(crate::EMAIL))
        .description(
            "
rosella recover is a tool for recovering MAGs from contigs using UMAP and HDBSCAN clustering.
            "
        );
    
    manual = manual.custom(binning_params_section());
    manual = manual.custom(reference_options_simple());
    manual = manual.custom(threads_options());

    manual = add_mapping_options(manual);
    manual = manual.custom(read_mapping_params_section());
    manual = add_thresholding_options(manual);
    manual = add_verbosity_flags(manual);
    
    manual
}

pub fn refine_full_help() -> Manual {
    let mut manual = Manual::new("rosella refine")
        .about(
            &format!(
                "Refine MAGs from contigs using UMAP and HDBSCAN clustering. (version {})",
                crate_version!()
            )
        )
        .author(Author::new(crate::AUTHOR).email(crate::EMAIL))
        .description(
            "
rosella refine is a tool for recovering MAGs from contigs using UMAP and HDBSCAN clustering.
            "
        );
    
    manual = manual.custom(binning_params_section());
    manual = manual.custom(reference_options_simple());
    manual = manual.custom(threads_options());
    manual = manual.custom(refining_options());

    manual = add_mapping_options(manual);
    manual = manual.custom(read_mapping_params_section());
    manual = add_thresholding_options(manual);
    manual = add_verbosity_flags(manual);
    
    manual
}
