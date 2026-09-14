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

    /// Pieces the protein file is cut into, each searched by its own hmmsearch. Defaults to a
    /// quarter of the thread count, since a shard costs a master thread plus its workers
    #[arg(long = "hmm-shards", value_parser = clap::value_parser!(u16).range(1..=64), hide_short_help = true)]
    pub hmm_shards: Option<u16>,

    /// Reuse the single copy marker annotation across runs over the same assembly, keyed on
    /// the build and every setting that changes it
    #[arg(long = "marker-cache", hide_short_help = true)]
    pub marker_cache: Option<String>,

    /// Write every single copy marker hit, with whether its gene ran off a contig end
    #[arg(long = "marker-report", hide_short_help = true)]
    pub marker_report: Option<String>,

    /// Write every candidate the rescue pool judged, with its rank, verdict and members
    #[arg(long = "pool-report", hide_short_help = true)]
    pub pool_report: Option<String>,

    /// Completeness a candidate needs before the pool adopts it
    #[arg(long = "min-completeness", default_value_t = crate::refine::rung::DEFAULT_COMPLETENESS, value_parser = crate::cli::common::percentage, hide_short_help = true)]
    pub min_completeness: f64,

    /// Contamination a candidate may carry before the pool refuses it
    #[arg(long = "max-contamination", default_value_t = crate::refine::rung::DEFAULT_CONTAMINATION, value_parser = crate::cli::common::percentage, hide_short_help = true)]
    pub max_contamination: f64,

    /// Whether every bin short of the bars goes back in the pot with the unbinned and is
    /// embedded again as one pool
    #[arg(long = "dissolve", value_parser = crate::recover::settings::DISSOLVE_NAMES, default_value = "on")]
    pub dissolve: String,

    /// Keep the refined bins as they are rather than offering the scorer whole bin pairs to
    /// fuse
    #[arg(long = "no-join", action = clap::ArgAction::SetTrue)]
    pub no_join: bool,

    /// Searches of the pool, each one over the whole of it, with the neighbour count halving
    /// each round so a genome the dense graph buries can still form its own community
    #[arg(long = "dissolve-rounds", default_value_t = 6, value_parser = clap::value_parser!(u16).range(1..=32), hide_short_help = true)]
    pub dissolve_rounds: u16,

    /// Cap on the passes over the pool, each one re-embedding what the pass before it left
    /// unclaimed. The passes stop on their own once one finds bins the model scores worse
    /// than the last
    #[arg(long = "dissolve-passes", default_value_t = 3, value_parser = clap::value_parser!(u16).range(1..=32), hide_short_help = true)]
    pub dissolve_passes: u16,

    /// Relax the accept bar one rung a pass, re-embedding between rungs, rather than walking
    /// every rung against one set of proposals
    #[arg(long = "dissolve-rung-per-pass", action = clap::ArgAction::SetTrue, hide_short_help = true)]
    pub dissolve_rung_per_pass: bool,

    /// Spend every pass over the pool rather than stopping at the first pass whose bins score
    /// worse than the last
    #[arg(long = "dissolve-all-passes", action = clap::ArgAction::SetTrue, hide_short_help = true)]
    pub dissolve_all_passes: bool,

    /// Share of the run's measured genome scale the pool's size floor sits at
    #[arg(long = "genome-floor-share", default_value_t = 1.0, value_parser = crate::cli::common::unit_interval, hide_short_help = true)]
    pub genome_floor_share: f64,

    /// Contig to genome map in CAMI binning format, offered to the pool as extra candidates.
    /// A probe: it asks whether the bar would take the right grouping if it were handed one
    #[arg(long = "dissolve-oracle", hide_short_help = true)]
    pub dissolve_oracle: Option<String>,

    /// Weight on contamination when ranking rescue candidates by worth
    #[arg(long = "worth-contamination", default_value_t = crate::refine::rung::DEFAULT_WORTH_CONTAMINATION, hide_short_help = true)]
    pub worth_contamination: f64,

    /// Contamination a bin may carry before worth charges it any
    #[arg(long = "worth-allowance", default_value_t = 0.0, value_parser = crate::cli::common::percentage, hide_short_help = true)]
    pub worth_allowance: f64,

    /// Completeness points below the bar a bin may sit and still draw contigs from the bins
    /// under --min-bin-size, which are discarded anyway
    #[arg(long = "recruit-near-bar", value_parser = crate::cli::common::percentage, hide_short_help = true)]
    pub recruit_near_bar: Option<f64>,

    /// Break equal worth on contig count, so a whole genome outranks a piece of it
    #[arg(long = "combine-size-tie", action = clap::ArgAction::SetTrue, hide_short_help = true)]
    pub combine_size_tie: bool,

    /// Partition seeds the ladder is built at. Every labelling from every seed reaches the
    /// per-bin combination
    #[arg(long = "partition-seeds", default_value_t = 3, value_parser = clap::value_parser!(u16).range(1..=16), hide_short_help = true)]
    pub partition_seeds: u16,

    /// Assembly graph in GFA format. Its links join the neighbour graph as extra edges
    #[arg(long = "assembly-graph")]
    pub assembly_graph: Option<String>,

    /// Weight an assembly graph link carries in the neighbour graph
    #[arg(long = "assembly-graph-weight", default_value_t = 0.75, value_parser = crate::cli::common::non_negative, hide_short_help = true)]
    pub assembly_graph_weight: f64,

    /// Completeness bar of the pool's last rung, as a share of the full bar
    #[arg(long = "rung-floor", default_value_t = crate::refine::rung::DEFAULT_RUNG_FLOOR, value_parser = crate::cli::common::unit_interval, hide_short_help = true)]
    pub rung_floor: f64,

    /// Ceiling on the contamination bar the pool's late rungs may relax to, as a multiple of
    /// the full bar
    #[arg(long = "rung-ceiling", default_value_t = crate::refine::rung::DEFAULT_RUNG_CEILING, value_parser = crate::cli::common::rung_ceiling, hide_short_help = true)]
    pub rung_ceiling: f64,

    /// Contamination a fused bin pair may carry before the merge is refused. Defaults to
    /// --max-contamination
    #[arg(long = "join-contamination", value_parser = crate::cli::common::percentage, hide_short_help = true)]
    pub join_contamination: Option<f64>,

    #[command(flatten)]
    pub full_help: FullHelp,

    #[command(flatten)]
    pub logging: Logging,
}
