use clap::{ArgGroup, Args};

use super::reports::MarkerReport;
use super::runtime::{Genomes, HelpFlags, Logging, Runtime};

#[derive(Args, Debug, Clone)]
#[command(disable_help_flag = true)]
#[command(group(ArgGroup::new("genomes").required(true).multiple(true)
    .args(["genome_fasta_files", "genome_fasta_directory"])))]
pub struct ScoreArgs {
    #[command(flatten)]
    pub genomes: Genomes,

    /// Where the quality table is written
    #[arg(short = 'o', long = "output-file", help_heading = "Input and output",
          default_value = crate::defaults::QUALITY_FILE)]
    pub output_file: String,

    #[command(flatten)]
    pub marker: MarkerReport,

    #[command(flatten)]
    pub runtime: Runtime,

    #[command(flatten)]
    pub logging: Logging,

    #[command(flatten)]
    pub help: HelpFlags,
}
