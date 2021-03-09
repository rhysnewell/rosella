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

const ALIGNMENT_OPTIONS: &'static str = "Define mapping(s) (required):
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
   -t, --threads <INT>                   Number of threads for mapping / sorting
                                         [default 8]
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

pub fn binning_full_help() -> &'static str {
    lazy_static! {
        static ref BINNING_HELP: String = format!(
    "rosella bin: Bins contigs from metagenomes into MAGs using coverage and TNF information

{}
{}

Binning parameters:
   -i, --coverage-values                 The output from the results of CoverM contig in MetaBAT mode
                                         on the provided assembly and short read samples. If not
                                         provided, rosella will calculate coverage values.
   -li, --longread-coverage-values       The output from the results of CoverM contig in MetaBAT mode
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
   -w, --window-size <FLOAT>             Window size in basepairs at which to calculate SNP and
                                         SV density. [default: 1000]
   --n-components <INT>                  Number of components for the UMAP algorithm to embed into. [default: 2]
   -n, --n-neighbors <INT>               Number of neighbors used in the UMAP algorithm. [default: 100]
   --a-spread <FLOAT>                    The spread of UMAP embeddings. Directly manipulates the
                                         \"a\" parameter. [default: 1.58]
   --b-tail <FLOAT>                      Similar to the heavy-tail parameter sometimes used in t-SNE.
                                         Directly manipulates the \"b\" parameter. [default: 0.4]
   --min-dist <FLOAT>                    Minimum dist parameter passed to UMAP algorithm. [default: 0.0]
   --scaler <STRING>                     Scaling method to use for coverage values and kmer frequencies.
                                         Options:
                                             - clr [default]
                                             - minmax
                                             - none


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
                                         log messages

Rhys J. P. Newell <rhys.newell near hdr.qut.edu.au>", ALIGNMENT_OPTIONS, MAPPER_HELP);
    }
    &BINNING_HELP
}

pub fn build_cli() -> App<'static, 'static> {
    // specify _2 lazily because need to define it at runtime.
    lazy_static! {


        static ref BINNING_HELP: String = format!(
            "
                            {}
              {}

{}

  rosella bin --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna --threads 10

{}

  coverm contig -m metabat --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna -o coverm.cov --threads 10
  rosella bin -i coverm.cov -r assembly.fna --output-directory rosella_out/ --threads 10

{}

  rosella bin -r assembly.fna --output-directory rosella_out/ --threads 10

See rosella bin --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "rosella bin"),
            ansi_term::Colour::Green.paint(
                "Recover MAGs from metagenomes using UMAP and HDBSCAN"),
            ansi_term::Colour::Purple.paint(
                "Example: Map paired reads to a reference to generate coverage info, then calculate kmer frequencies and bin MAGs"),
            ansi_term::Colour::Purple.paint(
                "Example: Use pregenerated coverage values, calculate kmer frequencies, and bin MAGs"),
            ansi_term::Colour::Purple.paint(
                "Example: Bin MAGs using pregenerated coverage and kmer results stored in output directory"),
        ).to_string();


    }

    return App::new("rosella")
        .version(crate_version!())
        .author("Rhys J.P. Newell <rhys.newell near hdr.qut.edu.au>")
        .about("MAG binner for metagenomes using coverage, TNF, and SNPs")
        .args_from_usage(
            "-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'",
        )
        .help(
            "
MAG binning using UMAP and HDBSCAN

Usage: rosella <subcommand> ...

Main subcommands:
\tbin \tMAG binning algorithm using coverage and TNF information across samples

Other options:
\t-V, --version\tPrint version information

Rhys J. P. Newell <r.newell near hdr.qut.edu.au>
",
        )
        .global_setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(
            SubCommand::with_name("bin")
                .about("Perform variant calling analysis and then binning")
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
                    Arg::with_name("variant-rates")
                        .long("variant-rates")
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
                        .default_value("./"),
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
                    Arg::with_name("parallel-contigs")
                        .long("parallel-contigs")
                        .default_value("8")
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
                    Arg::with_name("scaler")
                        .long("scaler")
                        .takes_value(true)
                        .possible_values(&["clr", "minmax", "none"])
                        .default_value("clr"),
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
                    Arg::with_name("coverage-fold")
                        .long("coverage-fold")
                        .default_value("1.0"),
                )
                .arg(
                    Arg::with_name("min-variant-depth")
                        .long("min-variant-depth")
                        .short("f")
                        .default_value("10"),
                )
                .arg(
                    Arg::with_name("min-variant-quality")
                        .long("min-variant-quality")
                        .default_value("10"),
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
                        .default_value("100"),
                )
                .arg(
                    Arg::with_name("n-components")
                        .long("n-components")
                        .default_value("2"),
                )
                .arg(
                    Arg::with_name("min-dist")
                        .long("min-dist")
                        .default_value("0.0"),
                )
                .arg(
                    Arg::with_name("a-spread")
                        .long("a-spread")
                        .default_value("1.58"),
                )
                .arg(Arg::with_name("b-tail").long("b-tail").default_value("0.4"))
                .arg(
                    Arg::with_name("minimum-reads-in-link")
                        .long("minimum-reads-in-link")
                        .default_value("5"),
                )
                .arg(
                    Arg::with_name("cluster-distance")
                        .long("cluster-distance")
                        .short("s")
                        .default_value("0.15"),
                )
                .arg(
                    Arg::with_name("base-quality-threshold")
                        .long("base-quality-threshold")
                        .short("q")
                        .default_value("13"),
                )
                .arg(
                    Arg::with_name("fdr-threshold")
                        .long("fdr-threshold")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("heterozygosity")
                        .long("heterozygosity")
                        .default_value("0.01"),
                )
                .arg(
                    Arg::with_name("indel-heterozygosity")
                        .long("indel-heterozygosity")
                        .default_value("0.001"),
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
                .arg(
                    Arg::with_name("window-size")
                        .long("window-size")
                        .short("w")
                        .default_value("1000"),
                )
                .arg(Arg::with_name("plot").long("plot"))
                .arg(Arg::with_name("include-longread-svs").long("include-longread-svs"))
                .arg(Arg::with_name("include-secondary").long("include-secondary"))
                .arg(Arg::with_name("include-soft-clipping").long("include-soft-clipping"))
                .arg(Arg::with_name("include-supplementary").long("include-supplementary"))
                .arg(
                    Arg::with_name("ploidy")
                        .long("ploidy")
                        .default_value("1")
                        .required(false),
                )
                .arg(
                    Arg::with_name("min-repeat-entropy")
                        .long("min-repeat-entropy")
                        .default_value("1.3")
                        .required(false),
                )
                .arg(Arg::with_name("force").long("force"))
                .arg(Arg::with_name("verbose").short("v").long("verbose"))
                .arg(Arg::with_name("quiet").long("quiet"))
                .arg(Arg::with_name("freebayes").long("freebayes"))
                .arg(
                    Arg::with_name("ulimit")
                        .long("ulimit")
                        .default_value("81920"),
                ),
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
