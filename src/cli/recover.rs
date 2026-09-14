use clap::{ArgAction, Args};

use super::binning::{BinningParams, DistanceParams, GraphParams, RefineParams};
use super::coverage::{
    AlignmentFlags, CoverageSource, CoverageTrimming, MappingParams, ReadFiltering,
};
use super::markers::MarkerParams;
use super::reports::ReportPaths;
use super::rescue::RescueParams;
use super::runtime::{Common, HelpFlags, Logging, Runtime, SeedParams, non_negative};

#[derive(Args, Debug, Clone)]
#[command(disable_help_flag = true)]
pub struct RecoverArgs {
    /// Assembly the contigs are read from
    #[arg(short = 'r', long, alias = "reference", help_heading = "Input and output")]
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

    /// Assembly graph in GFA format. Its links join the neighbour graph as extra edges
    #[arg(long = "assembly-graph", help_heading = "Neighbour graph")]
    pub assembly_graph: Option<String>,

    /// Weight an assembly graph link carries in the neighbour graph
    #[arg(long = "assembly-graph-weight", default_value_t = 0.75, value_parser = non_negative,
          requires = "assembly_graph", help_heading = "Neighbour graph", hide_short_help = true)]
    pub assembly_graph_weight: f64,

    #[command(flatten)]
    pub distance: DistanceParams,

    #[command(flatten)]
    pub refine: RefineParams,

    /// Keep the first clustering's bins instead of splitting the chimeric ones
    #[arg(long = "no-refine", action = ArgAction::SetTrue, help_heading = "Refinement")]
    pub no_refine: bool,

    /// Keep the refined bins as they are rather than offering the scorer whole bin pairs to
    /// fuse
    #[arg(long = "no-join", action = ArgAction::SetTrue, help_heading = "Refinement")]
    pub no_join: bool,

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
