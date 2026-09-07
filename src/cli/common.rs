use clap::{ArgAction, ArgGroup, Args};

use crate::clustering::clusterer::DEFAULT_LARGEST_CLUSTER;
use crate::clustering::graph_partition::{NODE_SIZE_NAMES, PARTITION_NAMES};
use crate::clustering::objective::OBJECTIVE_NAMES;
use crate::embedding::manifold::GRAPH_WEIGHT_NAMES;
use crate::embedding::metrics::{
    AGGREGATION_NAMES, BAND_NAMES, COMBINATION_NAMES, COMPOSITION_NAMES, VIEW_NAMES,
};
use crate::embedding::spectral::SPECTRAL_INIT_NAMES;
use crate::kmers::kmer_counting::DEFAULT_KMER_SIZE;
use crate::refine::bin_stats::SPLIT_LEVEL_NAMES;
use crate::refine::gates::SPLIT_GATE_NAMES;
use crate::refine::merger::MERGE_BAR_NAMES;
use crate::refine::solo::SOLO_POOL_NAMES;

/// Where coverage comes from. Any one of these is enough, so clap requires the group
/// rather than any single member.
#[derive(Args, Debug, Clone)]
#[command(group(ArgGroup::new("coverage-source").required(true).multiple(true).args([
    "read1", "coupled", "interleaved", "single", "longreads", "bam_files",
    "longread_bam_files", "coverage_file",
])))]
pub struct CoverageSource {
    /// Forward read files, paired with --read2
    #[arg(short = '1', long, num_args = 1.., action = ArgAction::Append, requires = "read2")]
    pub read1: Vec<String>,

    /// Reverse read files, paired with --read1
    #[arg(short = '2', long, num_args = 1.., action = ArgAction::Append, requires = "read1")]
    pub read2: Vec<String>,

    /// Paired read files, given as forward and reverse in turn
    #[arg(short = 'c', long, num_args = 1.., action = ArgAction::Append)]
    pub coupled: Vec<String>,

    /// Interleaved paired read files
    #[arg(long, num_args = 1.., action = ArgAction::Append)]
    pub interleaved: Vec<String>,

    /// Unpaired read files
    #[arg(long, num_args = 1.., action = ArgAction::Append)]
    pub single: Vec<String>,

    /// Long read files
    #[arg(long, num_args = 1.., action = ArgAction::Append)]
    pub longreads: Vec<String>,

    /// Reference sorted BAM files. No read mapping is undertaken
    #[arg(short = 'b', long = "bam-files", num_args = 1.., action = ArgAction::Append)]
    pub bam_files: Vec<String>,

    /// Reference sorted long read BAM files. No read mapping is undertaken
    #[arg(short = 'l', long = "longread-bam-files", num_args = 1.., action = ArgAction::Append)]
    pub longread_bam_files: Vec<String>,

    /// Precomputed CoverM coverage table, in place of mapping anything
    #[arg(short = 'C', long = "coverage-file")]
    pub coverage_file: Option<String>,
}

/// Passed straight to CoverM, which owns the list of mapper names and validates them.
#[derive(Args, Debug, Clone)]
pub struct MappingParams {
    /// Mapping software for short reads. Any name CoverM accepts
    #[arg(short = 'p', long)]
    pub mapper: Option<String>,

    /// Mapping software for long reads. Any name CoverM accepts
    #[arg(long = "longread-mapper")]
    pub longread_mapper: Option<String>,

    /// Extra parameters for minimap2, for indexing and for mapping. '-a' is always passed
    #[arg(
        long = "minimap2-parameters",
        alias = "minimap2-params",
        allow_hyphen_values = true
    )]
    pub minimap2_params: Option<String>,

    /// Extra parameters for BWA or BWA-MEM2
    #[arg(
        long = "bwa-parameters",
        alias = "bwa-params",
        allow_hyphen_values = true
    )]
    pub bwa_params: Option<String>,
}

#[derive(Args, Debug, Clone)]
pub struct ReadFiltering {
    /// Exclude reads aligned over fewer bases than this
    #[arg(long = "min-read-aligned-length")]
    pub min_read_aligned_length: Option<u32>,

    /// Exclude reads by percent identity to the reference, between 0 and 100
    #[arg(long = "min-read-percent-identity")]
    pub min_read_percent_identity: Option<f32>,

    /// Exclude reads by aligned percent of their length, between 0 and 100
    #[arg(long = "min-read-aligned-percent", default_value = "0.0")]
    pub min_read_aligned_percent: f32,

    /// As --min-read-aligned-length, but the pair is excluded together
    #[arg(long = "min-read-aligned-length-pair")]
    pub min_read_aligned_length_pair: Option<u32>,

    /// As --min-read-percent-identity, but the pair is excluded together
    #[arg(long = "min-read-percent-identity-pair")]
    pub min_read_percent_identity_pair: Option<f32>,

    /// As --min-read-aligned-percent, but the pair is excluded together
    #[arg(long = "min-read-aligned-percent-pair")]
    pub min_read_aligned_percent_pair: Option<f32>,
}

#[derive(Args, Debug, Clone)]
pub struct AlignmentFlags {
    /// Count only reads mapped in proper pairs
    #[arg(long = "proper-pairs-only", action = ArgAction::SetTrue)]
    pub proper_pairs_only: bool,

    /// Count secondary alignments
    #[arg(long = "include-secondary", action = ArgAction::SetTrue)]
    pub include_secondary: bool,

    /// Skip supplementary alignments
    #[arg(long = "exclude-supplementary", action = ArgAction::SetTrue)]
    pub exclude_supplementary: bool,
}

#[derive(Args, Debug, Clone)]
pub struct CoverageTrimming {
    /// Bases to ignore at each contig end
    #[arg(long = "contig-end-exclusion", default_value = "75")]
    pub contig_end_exclusion: usize,

    /// Discard this percent of the lowest coverage positions
    #[arg(long = "trim-min", default_value = "5.0")]
    pub trim_min: f32,

    /// Keep positions up to this percent of the coverage distribution
    #[arg(long = "trim-max", default_value = "95.0")]
    pub trim_max: f32,

    /// Report zero for contigs covered across less than this fraction
    #[arg(long = "min-covered-fraction", default_value = "0.0")]
    pub min_covered_fraction: f32,
}

#[derive(Args, Debug, Clone)]
pub struct BinningParams {
    /// Contigs shorter than this take no part in binning
    #[arg(long = "min-contig-size", default_value = "1500")]
    pub min_contig_size: usize,

    /// Clusters totalling less than this are not written as a bin
    #[arg(long = "min-bin-size", default_value = "200000")]
    pub min_bin_size: usize,

    /// Bins larger than this are always candidates for splitting
    #[arg(long = "max-bin-size", default_value = "15000000")]
    pub max_bin_size: usize,

    /// Neighbours per contig in the graph the embedding is built from
    #[arg(long = "n-neighbours", alias = "n-neighbors", default_value = "100")]
    pub n_neighbours: usize,

    /// Largest min_cluster_size the sweep tries
    #[arg(long = "max-cluster-size", value_parser = max_cluster_size_in_range,
          default_value_t = DEFAULT_LARGEST_CLUSTER)]
    pub max_cluster_size: usize,

    /// Rounds of refinement to attempt
    #[arg(long = "max-retries", default_value = "5")]
    pub max_retries: usize,

    /// What the parameter sweep ranks a labelling on
    #[arg(long = "objective", value_parser = OBJECTIVE_NAMES, default_value = "codelength")]
    pub objective: String,

    /// Where the cluster labels come from. The graph sources have no noise label, so every
    /// contig lands in a bin unless the eject takes it back out
    #[arg(long = "partition", value_parser = PARTITION_NAMES, default_value = "auto")]
    pub partition: String,

    /// What a node weighs in the partition: one per contig, or its length in bases with
    /// every edge scaled by the geometric mean of the two lengths it joins
    #[arg(long = "node-size", value_parser = NODE_SIZE_NAMES, default_value = "bp")]
    pub node_size: String,

    /// Pin the Leiden resolution instead of ranking a ladder of them on the objective
    #[arg(long = "partition-resolution", hide_short_help = true)]
    pub partition_resolution: Option<f64>,

    /// Sample the Leiden refinement target instead of taking the best gain. A fraction of
    /// the best gain available, so a larger value moves further from greedy
    #[arg(long = "partition-theta", value_parser = theta_above_zero, hide_short_help = true)]
    pub partition_theta: Option<f64>,

    /// Write every contig's nearest neighbours to this path and stop before embedding
    #[arg(long = "knn-report", hide_short_help = true)]
    pub knn_report: Option<std::path::PathBuf>,

    /// What a split has to clear. `validity` is the density validity alone, `floor` also wants
    /// two pieces at the bin floor, `genome` two at genome scale, `bimodal` two modes along the
    /// cut, and `auto` asks for the genome scale where the run can measure one and the modes
    /// where it cannot
    #[arg(long = "split-gate", value_parser = SPLIT_GATE_NAMES, default_value = "auto")]
    pub split_gate: String,

    /// Where the levels a bin is judged against come from. `derived` reads the run's own
    /// spread instead of flight's constants
    #[arg(long = "split-levels", value_parser = SPLIT_LEVEL_NAMES, default_value = "derived")]
    pub split_levels: String,

    /// Quantile of the run's own bin spreads a level sits at under `--split-levels derived`
    #[arg(long = "split-level-quantile", default_value = "0.75", value_parser = quantile_in_range)]
    pub split_level_quantile: f64,

    /// Also cut a bin in two on its own centroids, kept when the bin is bimodal along the cut
    #[arg(long = "bisect", action = clap::ArgAction::SetTrue)]
    pub bisect: bool,

    /// Keep contigs at least half the run's median closed genome together instead of standing
    /// each on its own when a bin holds more than one of them
    #[arg(long = "no-solo", action = clap::ArgAction::SetTrue)]
    pub no_solo: bool,

    /// Sweep every short contig into one leftover bin instead of offering each the genome-sized
    /// piece it sits nearest
    #[arg(long = "no-solo-scatter", action = clap::ArgAction::SetTrue)]
    pub no_solo_scatter: bool,

    /// Which contigs measure the run's genome scale. `alone` reads single-contig bins and
    /// unbinned contigs, `majority` also any contig holding over half its bin's bases, `long`
    /// every contig over the bin floor
    #[arg(long = "solo-pool", value_parser = SOLO_POOL_NAMES, default_value = "alone")]
    pub solo_pool: String,

    /// Keep a bin of one contig out of the merge. Without this it has no spread of its own
    /// to be judged by, so it is offered the spread its partner already tolerates
    #[arg(long = "no-merge-singles", action = clap::ArgAction::SetTrue)]
    pub no_merge_singles: bool,

    /// How loose a pair may be to merge. `pair` reads the two bins' own mean spread, `widest`
    /// the loosest contig each already holds, so the scale follows the contigs at hand
    #[arg(long = "merge-bar", value_parser = MERGE_BAR_NAMES, default_value = "pair")]
    pub merge_bar: String,

    /// Merge a pair only when each bin is the other's nearest, which asks for no distance at all
    #[arg(long = "merge-mutual", action = clap::ArgAction::SetTrue)]
    pub merge_mutual: bool,

    /// Merge only when one side holds less than the run's median closed genome
    #[arg(long = "merge-short-side", action = clap::ArgAction::SetTrue)]
    pub merge_short_side: bool,

    /// Compare contigs to each other with skani and keep two that align over most of both
    /// apart. Two loci of one genome do not align, the same locus in two organisms does
    #[arg(long = "homology", action = clap::ArgAction::SetTrue)]
    pub homology: bool,

    /// Cluster a bin again when two of its contigs align over most of both. Two organisms in
    /// one bin is a reason to re-cluster it, the way a duplicated single copy marker would be
    #[arg(long = "homology-trigger", action = clap::ArgAction::SetTrue)]
    pub homology_trigger: bool,

    /// Identity a pair has to reach to be called homologous
    #[arg(
        long = "homology-identity",
        default_value = "90.0",
        hide_short_help = true
    )]
    pub homology_identity: f64,

    /// Fraction of the shorter side's alignment a pair has to reach
    #[arg(
        long = "homology-aligned-fraction",
        default_value = "50.0",
        hide_short_help = true
    )]
    pub homology_aligned_fraction: f64,

    /// Length the longer contig of a pair has to reach, so two contigs short enough to be one
    /// repeat are not called two organisms. 0 asks nothing
    #[arg(
        long = "homology-min-length",
        default_value = "0",
        hide_short_help = true
    )]
    pub homology_min_length: usize,
}

fn theta_above_zero(value: &str) -> Result<f64, String> {
    let theta: f64 = value
        .parse()
        .map_err(|_| format!("`{value}` is not a number"))?;
    if theta > 0.0 && theta.is_finite() {
        Ok(theta)
    } else {
        Err(format!("`{value}` is not above 0"))
    }
}

fn knn_candidates_in_range(value: &str) -> Result<usize, String> {
    let candidates: usize = value
        .parse()
        .map_err(|_| format!("`{value}` is not a whole number"))?;
    if (2..=256).contains(&candidates) {
        Ok(candidates)
    } else {
        Err(format!("`{value}` is outside 2 to 256"))
    }
}

fn quantile_in_range(value: &str) -> Result<f64, String> {
    let quantile: f64 = value
        .parse()
        .map_err(|_| format!("`{value}` is not a number"))?;
    if (0.0..=1.0).contains(&quantile) {
        Ok(quantile)
    } else {
        Err(format!("`{value}` is outside [0, 1]"))
    }
}

/// Each carries a range because a typo would otherwise reach the optimiser as a useless
/// embedding.
#[derive(Args, Debug, Clone)]
pub struct EmbeddingOverrides {
    /// Dimensions in the embedding. Derived from the data's own dimensionality when unset
    #[arg(long = "n-components", value_parser = n_components_in_range)]
    pub n_components: Option<usize>,

    /// Neighbours and reverse neighbours each descent pass compares. Quadratic in the pass,
    /// so halving it quarters the work and loses recall
    #[arg(long = "knn-candidates", value_parser = knn_candidates_in_range, hide_short_help = true)]
    pub knn_candidates: Option<usize>,

    /// How the nearest neighbour graph becomes weighted edges. `fuzzy` is UMAP's smooth kNN
    /// sigma search, the other two read the neighbour lists directly and are far cheaper
    #[arg(long = "graph-weights", value_parser = GRAPH_WEIGHT_NAMES, default_value = "fuzzy",
          hide_short_help = true)]
    pub graph_weights: String,

    /// UMAP curve parameter a
    #[arg(long = "umap-a", value_parser = umap_a_in_range)]
    pub umap_a: Option<f32>,

    /// UMAP curve parameter b
    #[arg(long = "umap-b", value_parser = umap_b_in_range)]
    pub umap_b: Option<f32>,

    /// Smallest distance the layout packs points to. Fits the curve when set
    #[arg(long = "min-dist", value_parser = min_dist_in_range)]
    pub min_dist: Option<f32>,

    /// Scale of the embedded points. Fits the curve when set
    #[arg(long = "spread", value_parser = spread_in_range)]
    pub spread: Option<f32>,

    /// Layout optimisation epochs. Derived from the contig count when unset
    #[arg(long = "n-epochs", value_parser = n_epochs_in_range)]
    pub n_epochs: Option<usize>,

    /// Weight contig length into the graph edges
    #[arg(long = "length-weight", value_parser = length_weight_in_range, default_value = "0.0")]
    pub length_weight: f64,

    /// Which subspace the spectral start takes. `landmark` is seed free and scores far worse
    #[arg(long = "spectral-init", value_parser = SPECTRAL_INIT_NAMES,
          default_value = "random")]
    pub spectral_init: String,

    /// Report how much of each contig's neighbourhood the layout kept. Costs a second kNN
    /// build on every embedding, including one per bin during refinement
    #[arg(long = "report-preservation")]
    pub report_preservation: bool,
}

#[derive(Args, Debug, Clone)]
pub struct DistanceParams {
    /// How per-sample coverage distances combine
    #[arg(long = "coverage-aggregation", value_parser = AGGREGATION_NAMES,
          default_value = "arithmetic")]
    pub coverage_aggregation: String,

    /// Scale the variance floor by contig length
    #[arg(long = "length-scaled-variance", action = ArgAction::SetTrue)]
    pub length_scaled_variance: bool,

    /// Drop the coverage table's variance column and use the floor for every contig
    #[arg(long = "ignore-coverage-variance", action = ArgAction::SetTrue)]
    pub ignore_coverage_variance: bool,

    /// Views whose graphs are intersected. `combined` is the single distance
    #[arg(long = "embedding-views", value_parser = VIEW_NAMES, value_delimiter = ',',
          default_value = "combined")]
    pub embedding_views: Vec<String>,

    /// Coverage's share of the combined distance. Defaults to the samples that saw the pair,
    /// over that count plus one
    #[arg(long = "aggregate-weight", value_parser = unit_interval)]
    pub aggregate_weight: Option<f64>,

    /// Depth below this share of a pair's deepest sample counts as absent, and that sample is
    /// left out of the coverage distance
    #[arg(long = "presence-fraction", value_parser = unit_interval, default_value_t = 0.01)]
    pub presence_fraction: f64,

    /// How coverage and composition combine
    #[arg(long = "distance-combination", value_parser = COMBINATION_NAMES,
          default_value = "arithmetic")]
    pub distance_combination: String,

    /// Length of the k-mers the composition table counts
    #[arg(long = "kmer-size", value_parser = clap::value_parser!(u8).range(2..=6),
          default_value_t = DEFAULT_KMER_SIZE as u8)]
    pub kmer_size: u8,

    /// How two composition rows become one distance
    #[arg(long = "composition-metric", value_parser = COMPOSITION_NAMES,
          default_value = "rho")]
    pub composition_metric: String,

    /// Skip a sample in the coverage distance when its depth gap is inside the radius that
    /// already holds a neighbourhood. `drop` costs the sample its weight as well as its vote,
    /// `keep` costs only the vote
    #[arg(long = "coverage-band", value_parser = BAND_NAMES, default_value = "off")]
    pub coverage_band: String,
}

fn unit_interval(value: &str) -> Result<f64, String> {
    let parsed: f64 = value
        .parse()
        .map_err(|_| format!("`{value}` is not a number"))?;
    if (0.0..=1.0).contains(&parsed) {
        Ok(parsed)
    } else {
        Err(format!("`{parsed}` is outside 0.0 to 1.0"))
    }
}

/// Each stochastic stage draws from its own stream, so a run can hold three still and move
/// the fourth. Unset means the master seed, which is what keeps the default path unchanged.
#[derive(Args, Debug, Clone)]
pub struct SeedOverrides {
    /// Seed for the nearest neighbour graph. Defaults to --seed
    #[arg(long = "knn-seed", hide_short_help = true)]
    pub knn: Option<u64>,

    /// Seed for the spectral initialisation. Defaults to --seed
    #[arg(long = "init-seed", hide_short_help = true)]
    pub init: Option<u64>,

    /// Seed for the layout optimisation. Defaults to --seed
    #[arg(long = "layout-seed", hide_short_help = true)]
    pub layout: Option<u64>,

    /// Seed for the samples the objective and the refiner take. Defaults to --seed
    #[arg(long = "sample-seed", hide_short_help = true)]
    pub sample: Option<u64>,

    /// Seed for the node order a graph partition visits. Defaults to --seed
    #[arg(id = "partition-seed", long = "partition-seed", hide_short_help = true)]
    pub partition: Option<u64>,
}

#[derive(Args, Debug, Clone)]
pub struct Logging {
    /// Log at debug level
    #[arg(short, long, action = ArgAction::SetTrue)]
    pub verbose: bool,

    /// Log errors only
    #[arg(short, long, action = ArgAction::SetTrue)]
    pub quiet: bool,
}

#[derive(Args, Debug, Clone)]
pub struct Common {
    /// Where bins and the run's tables are written
    #[arg(short, long = "output-directory")]
    pub output_directory: String,

    /// Threads for the rayon pool and for CoverM
    #[arg(short, long, default_value = "10")]
    pub threads: usize,

    /// Seeds the embedding and every sample taken during clustering
    #[arg(long, default_value = "42")]
    pub seed: u64,

    /// Precomputed tetranucleotide frequency table, in place of counting them
    #[arg(short = 'K', long = "kmer-frequency-file")]
    pub kmer_frequency_file: Option<String>,
}

/// Rendered rather than parsed, so they are read off the command line before clap runs and
/// a missing required argument cannot stop the manual printing.
#[derive(Args, Debug, Clone)]
pub struct FullHelp {
    /// Print the full manual and exit
    #[arg(short = 'H', long = "full-help", action = ArgAction::SetTrue)]
    pub full_help: bool,

    /// Print the full manual as roff and exit
    #[arg(long = "full-help-roff", action = ArgAction::SetTrue)]
    pub full_help_roff: bool,
}

fn n_components_in_range(value: &str) -> Result<usize, String> {
    bounded(value, 2, 100)
}

fn n_epochs_in_range(value: &str) -> Result<usize, String> {
    bounded(value, 10, 10_000)
}

fn min_dist_in_range(value: &str) -> Result<f32, String> {
    bounded(value, 0.0, 5.0)
}

fn spread_in_range(value: &str) -> Result<f32, String> {
    bounded(value, 0.01, 10.0)
}

fn umap_a_in_range(value: &str) -> Result<f32, String> {
    bounded(value, 0.01, 10.0)
}

fn umap_b_in_range(value: &str) -> Result<f32, String> {
    bounded(value, 0.01, 5.0)
}

fn length_weight_in_range(value: &str) -> Result<f64, String> {
    bounded(value, 0.0, 2.0)
}

fn max_cluster_size_in_range(value: &str) -> Result<usize, String> {
    bounded(value, 2, 100_000)
}

pub(crate) fn eject_factor_in_range(value: &str) -> Result<f64, String> {
    bounded(value, 0.1, 10.0)
}

fn bounded<T>(value: &str, low: T, high: T) -> Result<T, String>
where
    T: std::str::FromStr + PartialOrd + std::fmt::Display + Copy,
{
    let parsed = value
        .parse::<T>()
        .map_err(|_| format!("`{value}` is not a number"))?;
    if parsed < low || parsed > high {
        return Err(format!("`{value}` is outside {low} to {high}"));
    }
    Ok(parsed)
}
