use clap::Args;

use super::coverage::ReadArgs;
use super::runtime::{HelpFlags, Logging, Runtime};

#[derive(Args, Debug, Clone)]
#[command(disable_help_flag = true)]
pub struct MakeArgs {
    /// Assembly to map reads to. BAMs need none
    #[arg(
        short = 'r',
        long,
        short_alias = 'a',
        alias = "reference",
        help_heading = "Input and output"
    )]
    pub assembly: Option<String>,

    /// Depth directory to write: contigs.tsv and each sample's 100 bp windows and table. Run
    /// once per sample or for all at once. recover and refine take it with --depth-directory
    #[arg(short, long = "output-directory", help_heading = "Input and output")]
    pub output_directory: String,

    #[command(flatten)]
    pub reads: ReadArgs,

    #[command(flatten)]
    pub runtime: Runtime,

    #[command(flatten)]
    pub logging: Logging,

    #[command(flatten)]
    pub help: HelpFlags,
}
