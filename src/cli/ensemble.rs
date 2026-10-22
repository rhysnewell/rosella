use clap::Args;

use super::recover::RecoverArgs;

#[derive(Args, Debug, Clone)]
#[command(disable_help_flag = true)]
pub struct EnsembleArgs {
    #[command(flatten)]
    pub recover: RecoverArgs,

    /// Bins, one entry per binner, gzipped or not. A directory holds one .fa, .fna or .fasta
    /// file per bin. A file is a table of contig and bin, tab separated, as a CAMI binning is
    #[arg(
        short = 'd',
        long = "bin-directories",
        num_args = 1..,
        required = true,
        help_heading = "Input and output"
    )]
    pub bin_directories: Vec<String>,
}
