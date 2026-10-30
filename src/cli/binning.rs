use clap::Args;

use crate::cli::runtime::{non_negative, percentage};
use crate::clustering::graph_partition::Partition;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Length {
    Auto,
    Given(usize),
}

impl Length {
    pub fn or(self, auto: usize) -> usize {
        match self {
            Self::Auto => auto,
            Self::Given(length) => length,
        }
    }
}

impl std::str::FromStr for Length {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        match text {
            "auto" => Ok(Self::Auto),
            _ => text
                .parse()
                .map(Self::Given)
                .map_err(|_| format!("expected auto or a length in bp, not {text}")),
        }
    }
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Binning")]
pub struct BinningParams {
    /// Shortest contig the partition takes. auto is 1500. Shorter contigs can still join a bin
    /// through --attach-floor
    #[arg(long = "min-contig-size", default_value = "auto")]
    pub min_contig_size: Length,

    /// Shortest contig that may join a bin after refinement, by where its neighbours sit. auto
    /// walks down from --min-contig-size. It stops where the joining contigs' markers say most
    /// of them are foreign, or where they would add less than a bin. A length given attaches
    /// from there at a neighbour share of one half
    #[arg(long = "attach-floor", default_value = "auto")]
    pub attach_floor: Length,

    /// Smallest bin written, in bp. auto is the 10 kb a replicon needs to be written alone. A
    /// contig on its own still needs the genome scale the run's near complete bins measure
    #[arg(long = "min-bin-size", default_value = "auto")]
    pub min_bin_size: Length,

    /// Bins over this many bp are always split candidates, unless they hold too few contigs to
    /// cluster again
    #[arg(long = "max-bin-size", default_value = "15000000")]
    pub max_bin_size: usize,

    #[arg(long = "partition", value_enum, default_value_t = Partition::Both,
          help = "Graph partition method. Neither leaves contigs out as noise. Under labelprop \
                  the rescue pool still runs Leiden, since it needs a ladder of rungs")]
    pub partition: Partition,

    /// Most training passes over the contigs for the learned space
    #[arg(long = "epochs", default_value_t = crate::defaults::EPOCHS,
          value_parser = clap::value_parser!(u32).range(1..))]
    pub epochs: u32,

    /// Annealed epochs without a real fall in loss before training stops. 0 trains for every
    /// --epochs
    #[arg(long = "patience", default_value_t = crate::defaults::PATIENCE)]
    pub patience: u32,

    /// Windows of each contig the learned space trains on beside the whole contig, up to 7.
    /// More windows cost training time. Read only with window depth
    #[arg(long = "stretches", default_value_t = crate::defaults::STRETCHES,
          value_parser = clap::value_parser!(u32).range(1..crate::learned::pieces::VIEWS as i64))]
    pub stretches: u32,

    /// Train the learned space mostly on contigs it has not learned yet. About a quarter of the
    /// training time, for a few bins fewer
    #[arg(long = "fast-train")]
    pub fast_train: bool,

    /// Train the learned space on extra steps whose contigs each carry a whole copy of one
    /// single copy marker, so are known to be different genomes. A carrier takes part when at
    /// least this percent of its reads tie with another contig. Read only with window depth
    #[arg(long = "marker-steps", value_parser = percentage)]
    pub marker_steps: Option<f64>,

    /// Train, embed and run the big neighbour searches on the GPU if there is one, else on the CPU
    #[arg(long = "gpu")]
    pub gpu: bool,

    /// Directory to keep the trained learned space in. A later recover or refine over the same
    /// assembly, depth and training flags reads it in place of training. Read only with window
    /// depth
    #[arg(long = "space-cache", hide_short_help = true)]
    pub space_cache: Option<String>,
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Neighbour graph")]
pub struct GraphParams {
    /// Neighbours per contig in the neighbour graph
    #[arg(long = "n-neighbours", alias = "n-neighbors", default_value = "100")]
    pub n_neighbours: usize,

    /// Neighbours of neighbours the search tries each round. Defaults to half of
    /// --n-neighbours. Too few and the seed decides the graph
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

impl BinningParams {
    pub fn cutoff(&self) -> usize {
        self.min_contig_size.or(crate::defaults::MIN_CONTIG_SIZE)
    }

    pub fn min_bin_size(&self) -> usize {
        self.min_bin_size.or(crate::markers::replicon::LEAST_BASES)
    }

    pub fn attach_floor(&self) -> Option<usize> {
        match self.attach_floor {
            Length::Auto => None,
            Length::Given(length) => Some(length.min(self.cutoff())),
        }
    }
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
    /// k-mer length for the composition table
    #[arg(long = "kmer-size", default_value_t = crate::kmers::kmer_counting::DEFAULT_KMER_SIZE,
          value_parser = clap::builder::RangedI64ValueParser::<usize>::new()
              .range(crate::kmers::kmer_counting::KMER_SIZES))]
    pub kmer_size: usize,

    /// Keep the composition table beside the bins for a later run over the same assembly. It
    /// runs to hundreds of MB on a large assembly and takes seconds to rebuild
    #[arg(long = "write-kmer-table", action = clap::ArgAction::SetTrue,
          hide_short_help = true)]
    pub write_kmer_table: bool,
}
