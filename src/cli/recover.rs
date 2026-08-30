use clap::Args;

use super::common::*;

#[derive(Args, Debug, Clone)]
pub struct RecoverArgs {
    /// Assembly the contigs are read from
    #[arg(short = 'r', long, alias = "reference")]
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
    pub overrides: EmbeddingOverrides,

    #[command(flatten)]
    pub distance: DistanceParams,

    #[command(flatten)]
    pub seeds: SeedOverrides,

    /// Split chimeric bins after the first clustering
    #[arg(long, action = clap::ArgAction::SetTrue)]
    pub refine: bool,

    /// Rejoin bins that are pieces of one genome
    #[arg(long, action = clap::ArgAction::SetTrue)]
    pub merge: bool,

    /// Re-cluster outliers on their own instead of offering them to the existing bins
    #[arg(long = "no-recruit", action = clap::ArgAction::SetTrue)]
    pub no_recruit: bool,

    #[command(flatten)]
    pub full_help: FullHelp,

    #[command(flatten)]
    pub logging: Logging,
}
