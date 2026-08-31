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

    /// Keep the first clustering's bins instead of splitting the chimeric ones
    #[arg(long = "no-refine", action = clap::ArgAction::SetTrue)]
    pub no_refine: bool,

    /// Keep bins that are pieces of one genome apart instead of rejoining them
    #[arg(long = "no-merge", action = clap::ArgAction::SetTrue)]
    pub no_merge: bool,

    /// Re-cluster outliers on their own instead of offering them to the existing bins
    #[arg(long = "no-recruit", action = clap::ArgAction::SetTrue)]
    pub no_recruit: bool,

    /// Keep contigs sitting further from their binmates than a typical bin's own spread
    #[arg(long = "no-eject", action = clap::ArgAction::SetTrue)]
    pub no_eject: bool,

    /// How far past the run's own level a contig has to sit before it is ejected
    #[arg(long = "eject-factor", default_value_t = 1.25, value_parser = crate::cli::common::eject_factor_in_range)]
    pub eject_factor: f64,

    #[command(flatten)]
    pub full_help: FullHelp,

    #[command(flatten)]
    pub logging: Logging,
}
