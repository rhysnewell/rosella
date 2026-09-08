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

    /// Keep the contigs a bin holds a second copy of, which is what two fused strains look like
    #[arg(long = "no-eject-duplicated", action = clap::ArgAction::SetTrue)]
    pub no_eject_duplicated: bool,

    /// Call single copy markers over the assembly and cut the bins that hold a second copy of
    /// a quarter of them on those markers. Measured at no gain on single sample: the markers
    /// name the fused bin but the graph holds no cut for the rest of its contigs to follow
    #[arg(long = "markers", action = clap::ArgAction::SetTrue)]
    pub markers: bool,

    /// Share of a bin's markers held twice before it is cut on them
    #[arg(long = "fusion-bar", default_value_t = crate::markers::DEFAULT_FUSION_BAR, value_parser = crate::cli::common::unit_interval, hide_short_help = true)]
    pub fusion_bar: f64,

    /// Keep the bins under the genome floor, and the ones holding their own sequence twice,
    /// where they are rather than embedding them again as one pool
    #[arg(long = "no-rescue", action = clap::ArgAction::SetTrue)]
    pub no_rescue: bool,

    /// Offer the contigs the rescue pool refused back to the bins that survived it
    #[arg(long = "recruit-rescued", action = clap::ArgAction::SetTrue)]
    pub recruit_rescued: bool,

    /// Duplicated share of a bin's k-mers before it is examined at all
    #[arg(long = "duplication-bar", default_value_t = crate::refine::duplication::DEFAULT_BAR, value_parser = crate::cli::common::unit_interval, hide_short_help = true)]
    pub duplication_bar: f64,

    /// Share of a contig's k-mers the rest of the bin must hold before it can leave
    #[arg(long = "duplication-link", default_value_t = crate::refine::duplication::DEFAULT_LINK, value_parser = crate::cli::common::unit_interval, hide_short_help = true)]
    pub duplication_link: f64,

    /// Sketch hashes a contig needs before its containment is trusted
    #[arg(long = "duplication-min-hashes", default_value_t = crate::refine::duplication::DEFAULT_MIN_HASHES, hide_short_help = true)]
    pub duplication_min_hashes: usize,

    /// k for the duplication sketch, which is not the composition k
    #[arg(long = "duplication-kmer-size", default_value_t = crate::kmers::sketch::DEFAULT_KMER_SIZE, value_parser = clap::value_parser!(u8).range(21..=31), hide_short_help = true)]
    pub duplication_kmer_size: u8,

    /// One k-mer in this many is kept in the sketch
    #[arg(long = "duplication-scale", default_value_t = crate::kmers::sketch::DEFAULT_SCALE, value_parser = clap::value_parser!(u64).range(1..=10000), hide_short_help = true)]
    pub duplication_scale: u64,

    #[command(flatten)]
    pub full_help: FullHelp,

    #[command(flatten)]
    pub logging: Logging,
}
