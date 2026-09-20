use clap::{ArgAction, Args};

use super::binning::{BinningParams, DistanceParams, GraphParams};
use super::coverage::{
    AlignmentFlags, CoverageSource, CoverageTrimming, MappingParams, ReadFiltering,
};
use super::markers::MarkerParams;
use super::reports::ReportPaths;
use super::rescue::RescueParams;
use super::runtime::{Common, HelpFlags, Logging, Runtime, SeedParams};

#[derive(Args, Debug, Clone)]
#[command(disable_help_flag = true)]
pub struct RecoverArgs {
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

    /// Keep the first clustering's bins instead of splitting the chimeric ones
    #[arg(long = "no-refine", action = ArgAction::SetTrue, help_heading = "Refinement")]
    pub no_refine: bool,

    /// Keep the refined bins as they are rather than offering the scorer whole bin pairs to
    /// fuse
    #[arg(long = "no-join", action = ArgAction::SetTrue, help_heading = "Refinement")]
    pub no_join: bool,

    /// Let a split through that leaves one genome standing and scatters less than a genome,
    /// rather than requiring it to leave two
    #[arg(long = "trim", action = ArgAction::SetTrue, help_heading = "Refinement")]
    pub trim: bool,

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
