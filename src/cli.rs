use clap::*;

const MAPPING_SOFTWARE_LIST: &[&str] = &[
    "bwa-mem",
    "minimap2-sr",
    "minimap2-ont",
    "minimap2-pb",
    "minimap2-no-preset",
    "ngmlr",
];
const DEFAULT_MAPPING_SOFTWARE: &str = "minimap2-sr";

const LONGREAD_MAPPING_SOFTWARE_LIST: &[&str] =
    &["minimap2-ont", "minimap2-pb", "ngmlr-ont", "ngmlr-pb"];
const DEFAULT_LONGREAD_MAPPING_SOFTWARE: &str = "minimap2-ont";

const MAPPER_HELP: &'static str = "  
    -p, --mapper <NAME>             Underlying mapping software used
                                    (\"minimap2-sr\", \"bwa-mem\",
                                    \"ngmlr-ont\", \"ngmlr-pb\", \"minimap2-ont\",
                                    \"minimap2-pb\", or \"minimap2-no-preset\").
                                    minimap2 -sr, -ont, -pb, -no-preset specify
                                    '-x' preset of minimap2 to be used
                                    (with map-ont, map-pb for -ont, -pb).
                                    [default: \"minimap2-sr\"] \n
    --minimap2-params PARAMS        Extra parameters to provide to minimap2,
                                    both indexing command (if used) and for
                                    mapping. Note that usage of this parameter
                                    has security implications if untrusted input
                                    is specified. '-a' is always specified.
                                    [default \"\"] \n
    --minimap2-reference-is-index   Treat reference as a minimap2 database, not
                                    as a FASTA file.\n
    --bwa-params PARAMS             Extra parameters to provide to BWA. Note
                                    that usage of this parameter has security
                                    implications if untrusted input is specified.
                                    [default \"\"]\n
    --ngmlr-params PARAMS           Extra parameters to provide to NGMLR.
                                    --bam-fix, -x ont, -t are already set. Note
                                    that usage of this parameter has security
                                    implications if untrusted input is specified.\n";

const ALIGNMENT_OPTIONS: &'static str = "Define mapping(s):
  Either define BAM:
   -b, --bam-files <PATH> ..             Path to BAM file(s). These must be
                                         reference sorted (e.g. with samtools sort)
                                         unless --sharded is specified, in which
                                         case they must be read name sorted (e.g.
                                         with samtools sort -n).
   -l, --longread-bam-files <PATH> ..    Path to BAM files(s) generated from longreads.
                                         Must be reference sorted.

  Or do mapping:
   -r, --reference <PATH> ..             FASTA file of contigs to be binned
   -1 <PATH> ..                          Forward FASTA/Q file(s) for mapping
   -2 <PATH> ..                          Reverse FASTA/Q file(s) for mapping
   -c, --coupled <PATH> <PATH> ..        One or more pairs of forward and reverse
                                         FASTA/Q files for mapping in order
                                         <sample1_R1.fq.gz> <sample1_R2.fq.gz>
                                         <sample2_R1.fq.gz> <sample2_R2.fq.gz> ..
   --interleaved <PATH> ..               Interleaved FASTA/Q files(s) for mapping.
   --single <PATH> ..                    Unpaired FASTA/Q files(s) for mapping.
   --longreads <PATH> ..                 pacbio or oxford nanopore long reads FASTA/Q files(s).
   -d, --bam-file-cache-directory        Directory to store cached BAM files. BAM files are stored
                                         in /tmp by default.";

const BINNING_OPTIONS: &'static str = "Binning parameters:
   -i, --coverage-values                 The output from the results of CoverM contig in MetaBAT mode
                                         on the provided assembly and short read samples. If not
                                         provided, rosella will calculate coverage values.
   --longread-coverage-values            The output from the results of CoverM contig in MetaBAT mode
                                         on the provided assembly and longread samples. If not provided, rosella
                                         will calculate coverage values.
   --kmer-frequencies                    The kmer frequency table created by rosella.
   --min-contig-size                     Minimum contig size in base pairs to be considered for binning.
                                         Contigs between 1000 bp and this value will be recovered in
                                         the contig rescue stage if multiple samples are available.
                                         [default: 1500]
   --min-bin-size                        Minimum bin size in base pairs for MAG to be reported. If a bin
                                         is smaller than this, then it will be split apart and the contigs
                                         will potentially be appended to another already established bin.
                                         [default: 200000]
   -k, --kmer-size <INT>                 K-mer size used to generate k-mer frequency
                                         table. [default: 4]
   -n, --n-neighbors <INT>               Number of neighbors used in the UMAP algorithm. [default: 100]


Alignment and contig filtering (optional):
   --min-read-aligned-length <INT>            Exclude reads with smaller numbers of
                                         aligned bases [default: 0]
   --min-read-percent-identity <FLOAT>        Exclude reads by overall percent
                                         identity e.g. 0.95 for 95%. [default 0.0]
   --min-read-aligned-percent <FLOAT>         Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned. [default 0.97]
   --min-read-aligned-length-pair <INT>       Exclude pairs with smaller numbers of
                                         aligned bases.
                                         Conflicts --proper-pairs-only. [default 0.0]
   --min-read-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95%.
                                         Conflicts --proper-pairs-only. [default 0.0]
   --min-read-aligned-percent-pair <FLOAT>    Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned.
                                         Conflicts --proper-pairs-only. [default 0.0]
   --proper-pairs-only                Allows reads to be mapped as improper pairs
   --include-supplementary               Includes read alignments flagged as supplementary
   --include-secondary                   Includes read alignments flagged as secondary



Other arguments (optional):
   -m, --method <METHOD>                 Method for calculating coverage.
                                         One or more (space separated) of:
                                           trimmed_mean
                                           mean
                                         A more thorough description of the different
                                         methods is available at
                                         https://github.com/rhysnewell/lorikeet
   -o, --output-directory <STRING>       Output directory for files. [default: output]
   --min-covered-fraction FRACTION       Contigs with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0.0]
   --contig-end-exclusion                Exclude bases at the ends of reference
                                         sequences from calculation [default: 75]
   --trim-min FRACTION                   Remove this smallest fraction of positions
                                         when calculating trimmed_mean
                                         [default: 0.05]
   --trim-max FRACTION                   Maximum fraction for trimmed_mean
                                         calculations [default: 0.95]
   -t, --threads                         Number of threads used. [default: 16]
   --no-zeros                            Omit printing of genomes that have zero
                                         coverage

   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages";

pub fn binning_full_help() -> &'static str {
    lazy_static! {
        static ref BINNING_HELP: String = format!(
    "rosella recover: Bins contigs from metagenomes into MAGs using coverage and TNF information
{}
{}
{}

Rhys J. P. Newell <rhys.newell near hdr.qut.edu.au>", BINNING_OPTIONS, ALIGNMENT_OPTIONS, MAPPER_HELP);
    }
    &BINNING_HELP
}

pub fn refining_full_help() -> &'static str {
    lazy_static! {
        static ref REFINING_HELP: String = format!(
    "rosella refine: Refine MAGs from previous binning attempts or other tools optionally guided by CheckM scores

Refining options:
   -a, --assembly <PATH>                 The original assembly used to produce the provided bins.
   -f, --genome-fasta-files <PATHS>      One or more paths to genome FASTA files to be refined
   -d, --genome-fasta-directory <PATH>   Directory containing FASTA files to be refined
   -x, --genome-fasta-extension <STR>    FASTA file extension in --genome-fasta-directory
                                         [default \"fna\"]
   --checkm-file                         CheckM1 or CheckM2 output table containing values for
                                         at least one of the input MAGs.
   --max-contamination                   Maximum valid bin contamination. Bins with contamination
                                         higher than this will be deconstructed. [default: 10]
   --contaminated-only                   Will limit the bin refining algorithm to only the MAGs
                                         that exceed the max contamination threshold
{}
{}
{}

Rhys J. P. Newell <rhys.newell near hdr.qut.edu.au>", BINNING_OPTIONS, ALIGNMENT_OPTIONS, MAPPER_HELP);
    }
    &REFINING_HELP
}

pub fn build_cli() -> App<'static, 'static> {
    // specify _2 lazily because need to define it at runtime.
    lazy_static! {


        static ref BINNING_HELP: String = format!(
            "
                            {}
              {}

{}

  rosella recover --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10

{}

  coverm contig -m metabat --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna -o coverm.cov --threads 10
  rosella recover -i coverm.cov -r assembly.fna --output-directory rosella_out/ --threads 10

{}

  rosella recover -r assembly.fna --output-directory rosella_out/ --threads 10

See rosella recover --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "rosella recover"),
            ansi_term::Colour::Green.paint(
                "Recover MAGs from metagenomes using UMAP and HDBSCAN"),
            ansi_term::Colour::Purple.paint(
                "Example: Map paired reads to a reference to generate coverage info, then calculate kmer frequencies and bin MAGs"),
            ansi_term::Colour::Purple.paint(
                "Example: Use pregenerated coverage values, calculate kmer frequencies, and bin MAGs"),
            ansi_term::Colour::Purple.paint(
                "Example: Bin MAGs using pregenerated coverage and kmer results stored in output directory"),
        ).to_string();


       static ref REFINING_HELP: String = format!(
            "
                            {}
              {}

{}

  rosella refine --assembly final_contigs.fasta -f contaminated_bin.1.fna --checkm-file checkm.out --threads 10 --coverage-values coverm.cov

{}

  rosella refine -a final_contigs.fasta -d contaminated_bins/ -x fna --checkm-file checkm.out --threads 10 --coverage-values coverm.cov

{}

  rosella refine -a final_contigs.fasta -d contaminated_bins/ -x fna --threads 10 --coverage-values coverm.cov

See rosella refine --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "rosella refine"),
            ansi_term::Colour::Green.paint(
                "Refine one or more MAGs via provided file paths and optional CheckM file"),
            ansi_term::Colour::Purple.paint(
                "Example: Refine a single MAG with CheckM results"),
            ansi_term::Colour::Purple.paint(
                "Example: Refine all bins present within directory with given file extension"),
            ansi_term::Colour::Purple.paint(
                "Example: The CheckM file is optional"),
        ).to_string();
    }

    return App::new("rosella")
        .version(crate_version!())
        .author("Rhys J.P. Newell <rhys.newell near hdr.qut.edu.au>")
        .about("MAG recovery algorithm for metagenomes using UMAP and HDBSCAN")
        .args_from_usage(
            "-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'",
        )
        .help(
            "
MAG recovery & refinement using UMAP and HDBSCAN

Usage: rosella <subcommand> ...

Main subcommands:
\trecover \tMAG recovery algorithm using coverage and TNF information across samples
\trefine \tRefine a given set of MAGs using the Rosella algorithm

Other options:
\t-V, --version\tPrint version information

Rhys J. P. Newell <r.newell near hdr.qut.edu.au>
",
        )
        .global_setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(
            SubCommand::with_name("recover")
                .about("Perform read mapping, coverage calculation, and binning")
                .help(BINNING_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(
                    Arg::with_name("coverage-values")
                        .short("i")
                        .long("coverage-values")
                        .takes_value(true)
                        .conflicts_with_all(&[
                            "read1",
                            "read2",
                            "coupled",
                            "interleaved",
                            "single",
                            "bam-files",
                            "full-help",
                        ]),
                )
                .arg(
                    Arg::with_name("longread-coverage-values")
                        // .short("c")
                        .long("longread-coverage-values")
                        .takes_value(true)
                        .conflicts_with_all(&["longreads", "longread-bam-files", "full-help"]),
                )
                .arg(
                    Arg::with_name("kmer-frequencies")
                        .long("kmer-frequencies")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("bam-files")
                        .short("b")
                        .long("bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .conflicts_with_all(&[
                            "read1",
                            "read2",
                            "coupled",
                            "interleaved",
                            "single",
                            "coverage-values",
                            "full-help",
                        ]),
                )
                .arg(
                    Arg::with_name("read1")
                        .short("-1")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read2")
                        .conflicts_with_all(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "coverage-values",
                            "single",
                            "full-help",
                        ]),
                )
                .arg(
                    Arg::with_name("read2")
                        .short("-2")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read1")
                        .conflicts_with_all(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "coverage-values",
                            "single",
                            "full-help",
                        ]),
                )
                .arg(
                    Arg::with_name("coupled")
                        .short("-c")
                        .long("coupled")
                        .multiple(true)
                        .takes_value(true)
                        .conflicts_with_all(&[
                            "bam-files",
                            "read1",
                            "interleaved",
                            "single",
                            "coverage-values",
                            "full-help",
                        ]),
                )
                .arg(
                    Arg::with_name("interleaved")
                        .long("interleaved")
                        .multiple(true)
                        .takes_value(true)
                        .conflicts_with_all(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "single",
                            "coverage-values",
                            "full-help",
                        ]),
                )
                .arg(
                    Arg::with_name("single")
                        .long("single")
                        .multiple(true)
                        .takes_value(true)
                        .conflicts_with_all(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "coverage-values",
                            "full-help",
                        ]),
                )
                .arg(
                    Arg::with_name("longreads")
                        .long("longreads")
                        .multiple(true)
                        .takes_value(true)
                        .required(false)
                        .conflicts_with_all(&["longread-bam-files"]),
                )
                .arg(
                    Arg::with_name("longread-bam-files")
                        .short("l")
                        .long("longread-bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required(false)
                        .conflicts_with_all(&["longreads"]),
                )
                .arg(
                    Arg::with_name("reference")
                        .short("r")
                        .long("reference")
                        .alias("genome-fasta-files")
                        .takes_value(true)
                        .required_unless_one(&["full-help"]),
                )
                .arg(
                    Arg::with_name("bam-file-cache-directory")
                        .long("bam-file-cache-directory")
                        .short("d")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("output-directory")
                        .long("output-directory")
                        .short("o")
                        .default_value("rosella_bins/"),
                )
                .arg(
                    Arg::with_name("vcfs")
                        .long("vcfs")
                        .multiple(true)
                        .required(false),
                )
                .arg(
                    Arg::with_name("threads")
                        .short("t")
                        .long("threads")
                        .default_value("16")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("mapper")
                        .short("p")
                        .long("mapper")
                        .possible_values(MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("longread-mapper")
                        .long("longread-mapper")
                        .possible_values(LONGREAD_MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_LONGREAD_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("minimap2-params")
                        .long("minimap2-params")
                        .long("minimap2-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index")
                        .requires("reference"),
                )
                .arg(
                    Arg::with_name("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                .arg(
                    Arg::with_name("discard-unmapped")
                        .long("discard-unmapped")
                        .requires("bam-file-cache-directory"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .default_value("0.0")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-contig-size")
                        .long("min-contig-size")
                        .takes_value(true)
                        .default_value("1500"),
                )
                .arg(
                    Arg::with_name("min-bin-size")
                        .long("min-bin-size")
                        .takes_value(true)
                        .default_value("200000"),
                )
                .arg(
                    Arg::with_name("min-samples")
                        .long("min-samples")
                        .takes_value(true)
                        .default_value("1"),
                )
                .arg(
                    Arg::with_name("method")
                        .short("m")
                        .long("method")
                        .takes_value(true)
                        .possible_values(&["trimmed_mean", "mean", "metabat"])
                        .default_value("trimmed_mean"),
                )
                .arg(
                    Arg::with_name("min-covered-fraction")
                        .long("min-covered-fraction")
                        .default_value("0.0"),
                )
                .arg(
                    Arg::with_name("mapq-threshold")
                        .long("mapq-threshold")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("n-neighbors")
                        .long("n-neighbors")
                        .short("n")
                        .default_value("200"),
                )
                .arg(
                    Arg::with_name("base-quality-threshold")
                        .long("base-quality-threshold")
                        .short("q")
                        .default_value("13"),
                )
                .arg(
                    Arg::with_name("kmer-size")
                        .long("kmer-size")
                        .short("k")
                        .default_value("4"),
                )
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("75"),
                )
                .arg(
                    Arg::with_name("trim-min")
                        .long("trim-min")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("trim-max")
                        .long("trim-max")
                        .default_value("0.95"),
                )
                .arg(Arg::with_name("no-zeros").long("no-zeros"))
                .arg(Arg::with_name("proper-pairs-only").long("proper-pairs-only"))
                .arg(Arg::with_name("plot").long("plot"))
                .arg(Arg::with_name("include-longread-svs").long("include-longread-svs"))
                .arg(Arg::with_name("include-secondary").long("include-secondary"))
                .arg(Arg::with_name("include-soft-clipping").long("include-soft-clipping"))
                .arg(Arg::with_name("include-supplementary").long("include-supplementary"))
                .arg(Arg::with_name("force").long("force"))
                .arg(Arg::with_name("verbose").short("v").long("verbose"))
                .arg(Arg::with_name("quiet").long("quiet")),
        )
        .subcommand(
            SubCommand::with_name("refine")
                .about("Refine a given set of bins using the Rosella algorithm")
                .help(REFINING_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(
                    Arg::with_name("coverage-values")
                        .short("i")
                        .long("coverage-values")
                        .takes_value(true)
                        .conflicts_with_all(&[
                            "read1",
                            "read2",
                            "coupled",
                            "interleaved",
                            "single",
                            "bam-files",
                            "full-help",
                        ]),
                )
                .arg(
                    Arg::with_name("longread-coverage-values")
                        // .short("c")
                        .long("longread-coverage-values")
                        .takes_value(true)
                        .conflicts_with_all(&["longreads", "longread-bam-files", "full-help"]),
                )
                .arg(
                    Arg::with_name("kmer-frequencies")
                        .long("kmer-frequencies")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("bam-files")
                        .short("b")
                        .long("bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .conflicts_with_all(&[
                            "read1",
                            "read2",
                            "coupled",
                            "interleaved",
                            "single",
                            "coverage-values",
                            "full-help",
                        ]),
                )
                .arg(
                    Arg::with_name("read1")
                        .short("-1")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read2")
                        .conflicts_with_all(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "coverage-values",
                            "single",
                            "full-help",
                        ]),
                )
                .arg(
                    Arg::with_name("read2")
                        .short("-2")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read1")
                        .conflicts_with_all(&[
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "coverage-values",
                            "single",
                            "full-help",
                        ]),
                )
                .arg(
                    Arg::with_name("coupled")
                        .short("-c")
                        .long("coupled")
                        .multiple(true)
                        .takes_value(true)
                        .conflicts_with_all(&[
                            "bam-files",
                            "read1",
                            "interleaved",
                            "single",
                            "coverage-values",
                            "full-help",
                        ]),
                )
                .arg(
                    Arg::with_name("interleaved")
                        .long("interleaved")
                        .multiple(true)
                        .takes_value(true)
                        .conflicts_with_all(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "single",
                            "coverage-values",
                            "full-help",
                        ]),
                )
                .arg(
                    Arg::with_name("single")
                        .long("single")
                        .multiple(true)
                        .takes_value(true)
                        .conflicts_with_all(&[
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "coverage-values",
                            "full-help",
                        ]),
                )
                .arg(
                    Arg::with_name("longreads")
                        .long("longreads")
                        .multiple(true)
                        .takes_value(true)
                        .required(false)
                        .conflicts_with_all(&["longread-bam-files"]),
                )
                .arg(
                    Arg::with_name("longread-bam-files")
                        .short("l")
                        .long("longread-bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required(false)
                        .conflicts_with_all(&["longreads"]),
                )
                .arg(
                    Arg::with_name("bam-file-cache-directory")
                        .long("bam-file-cache-directory")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("checkm-file")
                        .long("checkm-file")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("max-contamination")
                        .long("max-contamination")
                        .takes_value(true)
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("contaminated-only")
                        .long("contaminated-only")
                        .takes_value(false)
                        .requires("checkm-file")
                )
                .arg(
                    Arg::with_name("assembly")
                        .short("a")
                        .alias("reference")
                        .takes_value(true)
                        .required_unless_one(&["full-help"])
                )
                .arg(
                    Arg::with_name("genome-fasta-files")
                        .short("f")
                        .long("genome-fasta-files")
                        .multiple(true)
                        .conflicts_with("genome-fasta-directory")
                        .required_unless_one(&[
                            "genome-fasta-directory",
                            "full-help",
                        ])
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("genome-fasta-directory")
                        .short("d")
                        .long("genome-fasta-directory")
                        .conflicts_with("genome-fasta-files")
                        .required_unless_one(&[
                            "genome-fasta-files",
                            "full-help",
                        ])
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("genome-fasta-extension")
                        .long("genome-fasta-extension")
                        .short("x")
                        .default_value("fna")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("output-directory")
                        .long("output-directory")
                        .short("o")
                        .default_value("rosella_refined_bins/"),
                )
                .arg(
                    Arg::with_name("threads")
                        .short("t")
                        .long("threads")
                        .default_value("16")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("mapper")
                        .short("p")
                        .long("mapper")
                        .possible_values(MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("longread-mapper")
                        .long("longread-mapper")
                        .possible_values(LONGREAD_MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_LONGREAD_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("minimap2-params")
                        .long("minimap2-params")
                        .long("minimap2-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index")
                        .requires("reference"),
                )
                .arg(
                    Arg::with_name("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                .arg(
                    Arg::with_name("discard-unmapped")
                        .long("discard-unmapped")
                        .requires("bam-file-cache-directory"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .default_value("0.0")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .conflicts_with("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-contig-size")
                        .long("min-contig-size")
                        .takes_value(true)
                        .default_value("1500"),
                )
                .arg(
                    Arg::with_name("min-bin-size")
                        .long("min-bin-size")
                        .takes_value(true)
                        .default_value("0"),
                )
                .arg(
                    Arg::with_name("n-neighbors")
                        .long("n-neighbors")
                        .short("n")
                        .default_value("200"),
                )
                .arg(
                    Arg::with_name("method")
                        .short("m")
                        .long("method")
                        .takes_value(true)
                        .possible_values(&["trimmed_mean", "mean", "metabat"])
                        .default_value("trimmed_mean"),
                )
                .arg(
                    Arg::with_name("min-covered-fraction")
                        .long("min-covered-fraction")
                        .default_value("0.0"),
                )
                .arg(
                    Arg::with_name("mapq-threshold")
                        .long("mapq-threshold")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("kmer-size")
                        .long("kmer-size")
                        .short("k")
                        .default_value("4"),
                )
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("75"),
                )
                .arg(
                    Arg::with_name("trim-min")
                        .long("trim-min")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("trim-max")
                        .long("trim-max")
                        .default_value("0.95"),
                )
                .arg(Arg::with_name("no-zeros").long("no-zeros"))
                .arg(Arg::with_name("proper-pairs-only").long("proper-pairs-only"))
                .arg(Arg::with_name("include-secondary").long("include-secondary"))
                .arg(Arg::with_name("include-soft-clipping").long("include-soft-clipping"))
                .arg(Arg::with_name("include-supplementary").long("include-supplementary"))
                .arg(Arg::with_name("force").long("force"))
                .arg(Arg::with_name("verbose").short("v").long("verbose"))
                .arg(Arg::with_name("quiet").long("quiet")),
        )
        .subcommand(
            SubCommand::with_name("kmer")
                .about("Generate kmer count matrix for contigs")
                //                .help(CONTIG_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(
                    Arg::with_name("reference")
                        .short("-r")
                        .long("reference")
                        .takes_value(true)
                        .required(true),
                )
                .arg(Arg::with_name("verbose").short("v").long("verbose")),
        );
}
