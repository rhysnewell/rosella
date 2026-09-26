use clap::Args;

use crate::cli::runtime::non_negative;
use crate::clustering::graph_partition::PARTITION_NAMES;
use crate::clustering::leiden::NULL_NAMES;
use crate::kmers::kmer_counting::KmerSizes;

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Binning")]
pub struct BinningParams {
    /// Contigs shorter than this take no part in binning
    #[arg(long = "min-contig-size", default_value_t = crate::defaults::MIN_CONTIG_SIZE)]
    pub min_contig_size: usize,

    /// Clusters totalling less than this are not written as a bin
    #[arg(long = "min-bin-size", default_value = "200000")]
    pub min_bin_size: usize,

    /// Bins larger than this are split candidates whatever their spread, unless they hold
    /// too few contigs to re-cluster
    #[arg(long = "max-bin-size", default_value = "15000000")]
    pub max_bin_size: usize,

    /// Where the cluster labels come from. The graph sources have no noise label, so every
    /// contig lands in a bin unless the pool leaves it out. The rescue pool needs a ladder to
    /// walk its rungs over, so it runs Leiden even under labelprop
    #[arg(long = "partition", value_parser = PARTITION_NAMES, default_value = "both")]
    pub partition: String,

    /// Aim the resolution ladder at the bin size bounds rather than at fractions of the
    /// assembly's own mass
    /// Null model the partition scores against. cpm penalises community mass in bases,
    /// degree penalises it in edge weight
    #[arg(long = "leiden-null", value_parser = NULL_NAMES, default_value = "cpm", hide_short_help = true)]
    pub leiden_null: String,

    #[arg(long = "anchor-ladder", action = clap::ArgAction::SetTrue, hide_short_help = true)]
    pub anchor_ladder: bool,

    #[arg(long = "attractor-floor", hide_short_help = true)]
    pub attractor_floor: Option<usize>,
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Neighbour graph")]
pub struct GraphParams {
    /// Neighbours per contig in the graph the embedding is built from
    #[arg(long = "n-neighbours", alias = "n-neighbors", default_value = "100")]
    pub n_neighbours: usize,

    /// Neighbours of neighbours the descent tries each round. Defaults to half of
    /// --n-neighbours. Too few and the descent settles on a local optimum that the seed decides
    #[arg(long = "knn-candidates", hide_short_help = true)]
    pub knn_candidates: Option<usize>,

    /// Assembly graph in GFA format. Its links join the neighbour graph as extra edges
    #[arg(long = "assembly-graph")]
    pub assembly_graph: Option<String>,

    /// Weight an assembly graph link carries in the neighbour graph
    #[arg(long = "assembly-graph-weight", default_value_t = 0.75, value_parser = non_negative,
          requires = "assembly_graph", hide_short_help = true)]
    pub assembly_graph_weight: f64,
}

impl GraphParams {
    pub fn candidates(&self) -> usize {
        self.knn_candidates
            .unwrap_or_else(|| crate::embedding::knn::candidates(self.n_neighbours))
            .max(1)
    }
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Distances")]
pub struct DistanceParams {
    /// Length of the k-mers the composition table counts. A comma separated list concatenates
    /// one block per k, as in 2,3,4
    #[arg(long = "kmer-size", value_parser = clap::value_parser!(KmerSizes),
          default_value = "4")]
    pub kmer_size: KmerSizes,

    /// Keep the composition table beside the bins so a later run over the same assembly reuses
    /// it. It runs to hundreds of megabytes on a large assembly and costs seconds to rebuild
    #[arg(long = "write-kmer-table", action = clap::ArgAction::SetTrue,
          hide_short_help = true)]
    pub write_kmer_table: bool,

    /// Subtract the distance two contigs of their lengths would show anyway, so one threshold
    /// judges a short pair and a long pair alike
    #[arg(long = "calibrate-composition", action = clap::ArgAction::SetTrue,
          hide_short_help = true)]
    pub calibrate_composition: bool,

    #[arg(long = "weight-blocks", action = clap::ArgAction::SetTrue, hide_short_help = true)]
    pub weight_blocks: bool,
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Refinement")]
pub struct RefineParams {
    /// Rounds of refinement to attempt
    #[arg(long = "max-retries", default_value = "5")]
    pub max_retries: usize,
}
