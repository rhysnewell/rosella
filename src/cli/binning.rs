use clap::{ArgAction, Args};

use crate::cli::runtime::{above_zero, knn_candidates_in_range, unit_interval};
use crate::clustering::graph_partition::PARTITION_NAMES;
use crate::kmers::kmer_counting::{DEFAULT_KMER_SIZE, KMER_SIZES};

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Binning")]
pub struct BinningParams {
    /// Contigs shorter than this take no part in binning
    #[arg(long = "min-contig-size", default_value_t = crate::defaults::MIN_CONTIG_SIZE)]
    pub min_contig_size: usize,

    /// Clusters totalling less than this are not written as a bin
    #[arg(long = "min-bin-size", default_value = "200000")]
    pub min_bin_size: usize,

    /// Bins larger than this are always candidates for splitting
    #[arg(long = "max-bin-size", default_value = "15000000")]
    pub max_bin_size: usize,

    /// Where the cluster labels come from. The graph sources have no noise label, so every
    /// contig lands in a bin unless the pool leaves it out
    #[arg(long = "partition", value_parser = PARTITION_NAMES, default_value = "both")]
    pub partition: String,

    /// Pin the Leiden resolution instead of ranking a ladder of them on codelength
    #[arg(long = "partition-resolution", hide_short_help = true)]
    pub partition_resolution: Option<f64>,

    /// Sample the Leiden refinement target instead of taking the best gain. A fraction of
    /// the best gain available, so a larger value moves further from greedy
    #[arg(long = "partition-theta", value_parser = above_zero, hide_short_help = true)]
    pub partition_theta: Option<f64>,
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Neighbour graph")]
pub struct GraphParams {
    /// Neighbours per contig in the graph the embedding is built from
    #[arg(long = "n-neighbours", alias = "n-neighbors", default_value = "100")]
    pub n_neighbours: usize,

    /// Neighbours and reverse neighbours each descent pass compares. Quadratic in the pass,
    /// so halving it quarters the work and loses recall
    #[arg(long = "knn-candidates", value_parser = knn_candidates_in_range,
          default_value_t = crate::embedding::knn::MAX_CANDIDATES, hide_short_help = true)]
    pub knn_candidates: usize,
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Distances")]
pub struct DistanceParams {
    /// Length of the k-mers the composition table counts
    #[arg(long = "kmer-size", value_parser = clap::value_parser!(u8).range(KMER_SIZES),
          default_value_t = DEFAULT_KMER_SIZE as u8)]
    pub kmer_size: u8,

    /// Drop the coverage table's variance column and use the floor for every contig
    #[arg(long = "ignore-coverage-variance", action = ArgAction::SetTrue)]
    pub ignore_coverage_variance: bool,

    /// Depth below this share of a pair's deepest sample counts as absent, and that sample is
    /// left out of the coverage distance
    #[arg(long = "presence-fraction", value_parser = unit_interval, default_value_t = 0.01)]
    pub presence_fraction: f64,
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Refinement")]
pub struct RefineParams {
    /// Rounds of refinement to attempt
    #[arg(long = "max-retries", default_value = "5")]
    pub max_retries: usize,

    /// Quantile of the run's own bin spreads a level sits at
    #[arg(long = "split-level-quantile", default_value = "0.75", value_parser = unit_interval)]
    pub split_level_quantile: f64,

    /// Also cut a bin in two on its own centroids, kept when the bin is bimodal along the cut
    #[arg(long = "bisect", action = ArgAction::SetTrue)]
    pub bisect: bool,
}
