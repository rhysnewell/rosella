use clap::Args;

use super::recover::RecoverArgs;

#[derive(Args, Debug, Clone)]
#[command(disable_help_flag = true)]
pub struct EnsembleArgs {
    #[command(flatten)]
    pub recover: RecoverArgs,

    /// Directories of bins, one per binner. Every .fa, .fna or .fasta file in one is a bin,
    /// gzipped or not
    #[arg(
        short = 'd',
        long = "bin-directories",
        num_args = 1..,
        required = true,
        help_heading = "Input and output"
    )]
    pub bin_directories: Vec<String>,
}
