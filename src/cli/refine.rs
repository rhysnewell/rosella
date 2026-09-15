use clap::{ArgAction, ArgGroup, Args};

use super::binning::{BinningParams, DistanceParams, GraphParams, RefineParams};
use super::coverage::{
    AlignmentFlags, CoverageSource, CoverageTrimming, MappingParams, ReadFiltering,
};
use super::runtime::{Common, HelpFlags, Logging, Runtime, SeedParams};

#[derive(Args, Debug, Clone)]
#[command(disable_help_flag = true)]
#[command(group(ArgGroup::new("genomes").required(true).multiple(true)
    .args(["genome_fasta_files", "genome_fasta_directory"])))]
pub struct RefineArgs {
    /// Assembly the bins were built from
    #[arg(short = 'r', long, alias = "reference", help_heading = "Input and output")]
    pub assembly: String,

    /// Bins to refine
    #[arg(short = 'f', long = "genome-fasta-files", num_args = 1.., action = ArgAction::Append,
          help_heading = "Input and output")]
    pub genome_fasta_files: Vec<String>,

    /// Directory holding the bins to refine
    #[arg(short = 'd', long = "genome-fasta-directory", help_heading = "Input and output")]
    pub genome_fasta_directory: Option<String>,

    /// Extension of the bins inside --genome-fasta-directory
    #[arg(short = 'x', long = "genome-fasta-extension", help_heading = "Input and output",
          default_value = crate::defaults::FASTA_EXTENSION)]
    pub genome_fasta_extension: String,

    /// Written into the name of every bin this run produces
    #[arg(long = "bin-tag", default_value = "refined_1", help_heading = "Input and output")]
    pub bin_tag: String,

    #[command(flatten)]
    pub common: Common,

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
    pub graph: GraphParams,

    #[command(flatten)]
    pub distance: DistanceParams,

    #[command(flatten)]
    pub refine: RefineParams,

    /// Bin quality table, used to decide which bins to look at
    #[arg(long = "bin-quality", help_heading = "Refinement")]
    pub bin_quality: Option<String>,

    /// Bins over this contamination are always candidates for splitting
    #[arg(long = "split-contamination", default_value = "15.0", requires = "bin_quality",
          help_heading = "Refinement")]
    pub split_contamination: f64,

    /// Let a split through that leaves one genome standing and scatters less than a genome,
    /// rather than requiring it to leave two
    #[arg(long = "trim", action = clap::ArgAction::SetTrue, help_heading = "Refinement")]
    pub trim: bool,

    #[command(flatten)]
    pub seeds: SeedParams,

    #[command(flatten)]
    pub runtime: Runtime,

    #[command(flatten)]
    pub logging: Logging,

    #[command(flatten)]
    pub help: HelpFlags,
}
