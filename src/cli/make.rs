use clap::Args;

use super::coverage::{MappingParams, Samples};
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

    /// Depth directory to write: contigs.tsv, each sample's depth file and coverage.tsv.gz,
    /// the table of every sample in it. Run once per sample or for all at once, and with no
    /// reads to rebuild the table alone. recover and refine read it with --depth-directory
    #[arg(short, long = "output-directory", help_heading = "Input and output")]
    pub output_directory: String,

    #[command(flatten)]
    pub samples: Samples,

    #[command(flatten)]
    pub mapping: MappingParams,

    #[command(flatten)]
    pub runtime: Runtime,

    #[command(flatten)]
    pub logging: Logging,

    #[command(flatten)]
    pub help: HelpFlags,
}
