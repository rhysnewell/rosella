use clap::{ArgAction, ArgGroup, Args};

use crate::clustering::clusterer::DEFAULT_LARGEST_CLUSTER;
use crate::clustering::objective::OBJECTIVE_NAMES;
use crate::embedding::metrics::AGGREGATION_NAMES;

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
    #[arg(long = "objective", value_parser = OBJECTIVE_NAMES, default_value = "dbcv")]
    pub objective: String,
}

/// Pinning any of these bypasses the bounds the derived values are clamped to, so each
/// carries a range: a typo would otherwise reach the optimiser as a useless embedding.
#[derive(Args, Debug, Clone)]
pub struct EmbeddingOverrides {
    /// Dimensions in the embedding. Derived from the sample count when unset
    #[arg(long = "n-components", value_parser = n_components_in_range)]
    pub n_components: Option<usize>,

    /// UMAP curve parameter a. Derived when unset
    #[arg(long = "umap-a", value_parser = umap_a_in_range)]
    pub umap_a: Option<f32>,

    /// UMAP curve parameter b. Derived when unset
    #[arg(long = "umap-b", value_parser = umap_b_in_range)]
    pub umap_b: Option<f32>,

    /// Layout optimisation epochs. Derived from the contig count when unset
    #[arg(long = "n-epochs", value_parser = n_epochs_in_range)]
    pub n_epochs: Option<usize>,

    /// Weight contig length into the graph edges
    #[arg(long = "length-weight", value_parser = length_weight_in_range, default_value = "0.0")]
    pub length_weight: f64,
}

#[derive(Args, Debug, Clone)]
pub struct DistanceParams {
    /// How per-sample coverage distances combine
    #[arg(long = "coverage-aggregation", value_parser = AGGREGATION_NAMES,
          default_value = "geometric")]
    pub coverage_aggregation: String,

    /// Scale the variance floor by contig length
    #[arg(long = "length-scaled-variance", action = ArgAction::SetTrue)]
    pub length_scaled_variance: bool,
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
