use clap::Args;

use super::binning::{BinningParams, DistanceParams, GraphParams};
use super::coverage::ReadArgs;
use super::markers::MarkerParams;
use super::reports::{PartitionProbes, ReportPaths};
use super::rescue::RescueParams;
use super::runtime::{Common, HelpFlags, Logging, Runtime, SeedParams};

// Everything recover and refine both read, so the two subcommands cannot drift apart.
#[derive(Args, Debug, Clone)]
pub struct RunArgs {
    /// Assembly the contigs are read from
    #[arg(
        short = 'r',
        long,
        short_alias = 'a',
        alias = "reference",
        help_heading = "Input and output"
    )]
    pub assembly: String,

    #[command(flatten)]
    pub common: Common,

    #[command(flatten)]
    pub reads: ReadArgs,

    #[command(flatten)]
    pub binning: BinningParams,

    #[command(flatten)]
    pub graph: GraphParams,

    #[command(flatten)]
    pub distance: DistanceParams,

    #[command(flatten)]
    pub rescue: RescueParams,

    #[command(flatten)]
    pub markers: MarkerParams,

    #[command(flatten)]
    pub reports: ReportPaths,

    #[command(flatten)]
    pub seeds: SeedParams,

    #[command(flatten)]
    pub runtime: Runtime,

    #[command(flatten)]
    pub logging: Logging,

    #[command(flatten)]
    pub help: HelpFlags,
}

#[derive(Args, Debug, Clone)]
#[command(disable_help_flag = true)]
pub struct RecoverArgs {
    #[command(flatten)]
    pub run: RunArgs,

    /// Seeds the first partition's ladder is built at. Each seed and arm gives its best rung
    /// to the per-bin combination
    #[arg(long = "partition-seeds", default_value_t = 2, help_heading = "Binning",
          value_parser = clap::value_parser!(u16).range(1..=16), hide_short_help = true)]
    pub partition_seeds: u16,

    #[command(flatten)]
    pub probes: PartitionProbes,
}
