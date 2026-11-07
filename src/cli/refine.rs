use clap::{ArgGroup, Args};

use super::recover::RunArgs;
use super::runtime::Genomes;

#[derive(Args, Debug, Clone)]
#[command(disable_help_flag = true)]
#[command(group(ArgGroup::new("genomes").required(true).multiple(true)
    .args(["genome_fasta_files", "genome_fasta_directory"])))]
pub struct RefineArgs {
    #[command(flatten)]
    pub genomes: Genomes,

    /// Tag in the name of every bin this run writes
    #[arg(
        long = "bin-tag",
        default_value = "refined_1",
        help_heading = "Input and output"
    )]
    pub bin_tag: String,

    #[command(flatten)]
    pub run: RunArgs,

    /// Bin quality table from CheckM1, CheckM2 or AMBER. A bin over --split-contamination in
    /// it is always a split candidate
    #[arg(
        long = "bin-quality",
        alias = "checkm-file",
        help_heading = "Refinement"
    )]
    pub bin_quality: Option<String>,

    /// Contamination in --bin-quality over which a bin is always a split candidate
    #[arg(
        long = "split-contamination",
        default_value = "15.0",
        requires = "bin_quality",
        help_heading = "Refinement"
    )]
    pub split_contamination: f64,
}
