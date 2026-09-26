use clap::{ArgAction, ArgGroup, Args};

use super::markers::MarkerParams;
use super::runtime::{HelpFlags, Logging, Runtime};

#[derive(Args, Debug, Clone)]
#[command(disable_help_flag = true)]
#[command(group(ArgGroup::new("genomes").required(true).multiple(true)
    .args(["genome_fasta_files", "genome_fasta_directory"])))]
pub struct ScoreArgs {
    /// Assembly the bins were built from
    #[arg(
        short = 'r',
        long,
        alias = "reference",
        help_heading = "Input and output"
    )]
    pub assembly: String,

    /// Bins to score
    #[arg(short = 'f', long = "genome-fasta-files", num_args = 1.., action = ArgAction::Append,
          help_heading = "Input and output")]
    pub genome_fasta_files: Vec<String>,

    /// Directory holding the bins to score
    #[arg(
        short = 'd',
        long = "genome-fasta-directory",
        help_heading = "Input and output"
    )]
    pub genome_fasta_directory: Option<String>,

    /// Extension of the bins inside --genome-fasta-directory
    #[arg(short = 'x', long = "genome-fasta-extension", help_heading = "Input and output",
          default_value = crate::defaults::FASTA_EXTENSION)]
    pub genome_fasta_extension: String,

    /// Where the quality table is written
    #[arg(short = 'o', long = "output-file", help_heading = "Input and output",
          default_value = crate::defaults::QUALITY_FILE)]
    pub output_file: String,

    /// Contigs shorter than this are left out, as they are in a recover run
    #[arg(long = "min-contig-size", help_heading = "Binning",
          default_value_t = crate::defaults::MIN_CONTIG_SIZE)]
    pub min_contig_size: usize,

    #[command(flatten)]
    pub markers: MarkerParams,

    #[command(flatten)]
    pub runtime: Runtime,

    #[command(flatten)]
    pub logging: Logging,

    #[command(flatten)]
    pub help: HelpFlags,
}
