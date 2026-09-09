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

    /// Keep the contigs a bin holds a second copy of, which is what two fused strains look like
    #[arg(long = "no-eject-duplicated", action = clap::ArgAction::SetTrue)]
    pub no_eject_duplicated: bool,

    /// Protein database for the gene family search. Without it, and without CHECKM2DB set, the
    /// pool falls back to scoring a candidate on the sequence it holds twice
    #[arg(long = "checkm2-db")]
    pub checkm2_db: Option<String>,

    /// Directory to keep the gene family tables in, so a rerun over the same assembly and
    /// database skips the search
    #[arg(long = "checkm2-cache")]
    pub checkm2_cache: Option<String>,

    /// Completeness a candidate needs before the pool adopts it
    #[arg(long = "min-completeness", default_value_t = crate::refine::rung::DEFAULT_COMPLETENESS, value_parser = crate::cli::common::percentage, hide_short_help = true)]
    pub min_completeness: f64,

    /// Contamination a candidate may carry before the pool refuses it
    #[arg(long = "max-contamination", default_value_t = crate::refine::rung::DEFAULT_CONTAMINATION, value_parser = crate::cli::common::percentage, hide_short_help = true)]
    pub max_contamination: f64,

    /// Keep the bins under the genome floor, and the ones holding their own sequence twice,
    /// where they are rather than embedding them again as one pool
    #[arg(long = "no-dissolve", action = clap::ArgAction::SetTrue)]
    pub no_dissolve: bool,

    /// Searches of the pool, each one over the whole of it, with the neighbour count halving
    /// each round so a genome the dense graph buries can still form its own community
    #[arg(long = "dissolve-rounds", default_value_t = 6, value_parser = clap::value_parser!(u16).range(1..=32), hide_short_help = true)]
    pub dissolve_rounds: u16,

    /// Cap on the passes over the pool, each one re-embedding what the pass before it left
    /// unclaimed. The passes stop on their own once one finds bins the model scores worse
    /// than the last
    #[arg(long = "dissolve-passes", default_value_t = 3, value_parser = clap::value_parser!(u16).range(1..=32), hide_short_help = true)]
    pub dissolve_passes: u16,

    /// Reuse the pool's first neighbour build for its later passes instead of rebuilding.
    /// Around a fifth off the wall for roughly one bin in a hundred and seventy
    #[arg(long = "fast-pool", action = clap::ArgAction::SetTrue)]
    pub fast_pool: bool,

    /// Contig to genome map in CAMI binning format, offered to the pool as extra candidates.
    /// A probe: it asks whether the bar would take the right grouping if it were handed one
    #[arg(long = "dissolve-oracle", hide_short_help = true)]
    pub dissolve_oracle: Option<String>,

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
