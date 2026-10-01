use clap::{ArgAction, ArgGroup, Args};

use super::runtime::{HelpFlags, Logging, Runtime};

#[derive(Args, Debug, Clone)]
#[command(disable_help_flag = true)]
#[command(group(ArgGroup::new("genomes").required(true).multiple(true)
    .args(["genome_fasta_files", "genome_fasta_directory"])))]
pub struct ScoreArgs {
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

    /// Write every single copy marker hit, with whether its gene ran off a contig end
    #[arg(
        long = "marker-report",
        help_heading = "Reports",
        hide_short_help = true
    )]
    pub marker_report: Option<String>,

    #[command(flatten)]
    pub runtime: Runtime,

    #[command(flatten)]
    pub logging: Logging,

    #[command(flatten)]
    pub help: HelpFlags,
}
