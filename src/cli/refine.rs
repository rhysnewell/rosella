use clap::{ArgAction, ArgGroup, Args};

use super::common::*;

#[derive(Args, Debug, Clone)]
#[command(group(ArgGroup::new("genomes").required(true).multiple(true)
    .args(["genome_fasta_files", "genome_fasta_directory"])))]
pub struct RefineArgs {
    /// Assembly the bins were built from. Not needed when both tables are supplied
    #[arg(short = 'r', long, alias = "reference",
          required_unless_present_all = ["coverage_file", "kmer_frequency_file"])]
    pub assembly: Option<String>,

    #[command(flatten)]
    pub common: Common,

    /// Bins to refine
    #[arg(short = 'f', long = "genome-fasta-files", num_args = 1.., action = ArgAction::Append)]
    pub genome_fasta_files: Vec<String>,

    /// Directory holding the bins to refine
    #[arg(short = 'd', long = "genome-fasta-directory")]
    pub genome_fasta_directory: Option<String>,

    /// Extension of the bins inside --genome-fasta-directory
    #[arg(short = 'x', long = "genome-fasta-extension", default_value = "fna")]
    pub genome_fasta_extension: String,

    /// CheckM1, CheckM2 or AMBER table, used to decide which bins to look at
    #[arg(long = "checkm-results")]
    pub checkm_results: Option<String>,

    /// Bins over this contamination are always candidates for splitting
    #[arg(long = "max-contamination", default_value = "15.0")]
    pub max_contamination: f64,

    /// Bins with fewer contigs than this are passed through untouched
    #[arg(long = "min-contig-count", default_value = "10")]
    pub min_contig_count: usize,

    /// Written into the name of every bin this run produces
    #[arg(long = "bin-tag", default_value = "refined_1")]
    pub bin_tag: String,

    #[command(flatten)]
    pub coverage: CoverageSource,

    #[command(flatten)]
    pub mapping: MappingParams,

    #[command(flatten)]
    pub filtering: ReadFiltering,

    #[command(flatten)]
    pub alignment: AlignmentFlags,

    #[command(flatten)]
    pub trimming: CoverageTrimming,

    #[command(flatten)]
    pub binning: BinningParams,

    #[command(flatten)]
    pub overrides: EmbeddingOverrides,

    #[command(flatten)]
    pub distance: DistanceParams,

    #[command(flatten)]
    pub full_help: FullHelp,

    #[command(flatten)]
    pub logging: Logging,
}
