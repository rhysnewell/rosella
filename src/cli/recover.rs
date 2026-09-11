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

    /// Deprecated, due for removal. Judge candidates on gene families from this protein
    /// database rather than on the single copy markers built into the binary
    #[arg(long = "gene-database")]
    pub gene_database: Option<String>,

    /// Search a gene cut by a contig end again, at the share of the model it could still match
    #[arg(long = "marker-fragments", action = clap::ArgAction::SetTrue, hide_short_help = true)]
    pub marker_fragments: bool,

    /// Judge a marker bin on its markers alone, without the run's genome scale floor
    #[arg(long = "marker-no-scale-floor", action = clap::ArgAction::SetTrue, hide_short_help = true)]
    pub marker_no_scale_floor: bool,

    /// Least share of a model a cut gene must span before its hit is trusted
    #[arg(long = "marker-fragment-span", default_value_t = crate::markers::fragments::DEFAULT_SPAN, hide_short_help = true)]
    pub marker_fragment_span: f64,

    /// How much lower the marker bar sits than the requested completeness
    #[arg(long = "marker-bar-offset", default_value_t = crate::markers::DEFAULT_BAR_OFFSET, value_parser = crate::cli::common::percentage, hide_short_help = true)]
    pub marker_bar_offset: f64,

    /// Contigs shorter than this are not searched for genes, though they are still binned
    #[arg(long = "gene-min-length", default_value = "0", hide_short_help = true)]
    pub gene_min_length: usize,

    /// Run the full path search for only this many metagenomic models, ranked on their best
    /// node. 0 runs every model the GC window admits
    #[arg(long = "gene-model-depth", default_value = "0", hide_short_help = true)]
    pub gene_model_depth: usize,

    /// Pieces the protein file is cut into, each searched by its own hmmsearch. Defaults to
    /// half the thread count, since a shard costs a master thread plus at least one worker
    #[arg(long = "hmm-shards", value_parser = clap::value_parser!(u16).range(1..=64), hide_short_help = true)]
    pub hmm_shards: Option<u16>,

    /// Write every single copy marker hit, with whether its gene ran off a contig end
    #[arg(long = "marker-report", hide_short_help = true)]
    pub marker_report: Option<String>,

    /// Write every candidate the rescue pool judged, with its rank, verdict and members
    #[arg(long = "pool-report", hide_short_help = true)]
    pub pool_report: Option<String>,

    /// Deprecated, due for removal with --gene-database. Directory to keep the gene family
    /// tables in, so a rerun over the same assembly and database skips the search
    #[arg(long = "gene-cache")]
    pub gene_cache: Option<String>,

    /// Deprecated, due for removal with --gene-database. Search the assembly again rather
    /// than reading or writing the gene family tables
    #[arg(long = "no-gene-cache", action = clap::ArgAction::SetTrue)]
    pub no_gene_cache: bool,

    /// Deprecated, due for removal with --gene-database. How hard the gene family search
    /// looks. The faster tiers drop the weakest hits, which the model reads as absent genes
    #[arg(long = "gene-sensitivity", default_value = "default", value_parser = ["default", "fast", "faster"], hide_short_help = true)]
    pub gene_sensitivity: String,

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

    /// Keep the refined bins as they are rather than offering the scorer whole bin pairs to
    /// fuse
    #[arg(long = "no-join", action = clap::ArgAction::SetTrue)]
    pub no_join: bool,

    /// Build the rescue pool's candidates from the graph alone, without the merge order the
    /// markers propose
    #[arg(long = "no-linkage", action = clap::ArgAction::SetTrue)]
    pub no_linkage: bool,

    /// Searches of the pool, each one over the whole of it, with the neighbour count halving
    /// each round so a genome the dense graph buries can still form its own community
    #[arg(long = "dissolve-rounds", default_value_t = 6, value_parser = clap::value_parser!(u16).range(1..=32), hide_short_help = true)]
    pub dissolve_rounds: u16,

    /// Cap on the passes over the pool, each one re-embedding what the pass before it left
    /// unclaimed. The passes stop on their own once one finds bins the model scores worse
    /// than the last
    #[arg(long = "dissolve-passes", default_value_t = 3, value_parser = clap::value_parser!(u16).range(1..=32), hide_short_help = true)]
    pub dissolve_passes: u16,

    /// Rebuild the pool's neighbours for every pass rather than inducing the later
    /// passes from the first build, which costs about a fifth of the wall
    #[arg(long = "no-fast-pool", action = clap::ArgAction::SetTrue, hide_short_help = true)]
    pub no_fast_pool: bool,

    /// Contig to genome map in CAMI binning format, offered to the pool as extra candidates.
    /// A probe: it asks whether the bar would take the right grouping if it were handed one
    #[arg(long = "dissolve-oracle", hide_short_help = true)]
    pub dissolve_oracle: Option<String>,

    /// Duplicated share of a bin's k-mers before it is examined at all
    #[arg(long = "duplication-bar", default_value_t = crate::refine::rung::DEFAULT_DUPLICATION_BAR, value_parser = crate::cli::common::unit_interval, hide_short_help = true)]
    pub duplication_bar: f64,


    /// Sketch hashes a contig needs before its containment is trusted
    /// Weight on contamination when ranking rescue candidates by worth
    #[arg(long = "worth-contamination", default_value_t = crate::refine::rung::DEFAULT_WORTH_CONTAMINATION, hide_short_help = true)]
    pub worth_contamination: f64,

    /// Completeness bar of the pool's last rung, as a share of the full bar
    #[arg(long = "rung-floor", default_value_t = crate::refine::rung::DEFAULT_RUNG_FLOOR, value_parser = crate::cli::common::unit_interval, hide_short_help = true)]
    pub rung_floor: f64,



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
