use std::path::PathBuf;

use clap::Args;

// None of these changes a bin. They write what a run decided, or hand it something to judge.
#[derive(Args, Debug, Clone, Default)]
#[command(next_help_heading = "Reports")]
pub struct ReportPaths {
    #[command(flatten)]
    pub marker: MarkerReport,

    /// Write every candidate the rescue pool judged, with its rank, verdict and members
    #[arg(long = "pool-report", hide_short_help = true)]
    pub pool_report: Option<PathBuf>,

    /// Write every parked short contig with the bin attach would put it in, its neighbour
    /// share, whether its markers repeat there and whether the bin refused it
    #[arg(long = "attach-report", hide_short_help = true)]
    pub attach_report: Option<PathBuf>,

    /// Write every contig the shed would unbin against every bin that holds a neighbour of
    /// it, with the rival's share, the claim and whether the contig completes it
    #[arg(long = "shed-report", hide_short_help = true)]
    pub shed_report: Option<PathBuf>,

    /// Write every contig a stage takes out of a bin, with the final bin its old bin became and
    /// what handing it back would do to that bin's completeness and contamination
    #[arg(long = "cut-report", hide_short_help = true)]
    pub cut_report: Option<PathBuf>,

    #[arg(
        long = "points-report",
        hide_short_help = true,
        help = "Write every contig read before the partition with its point in the learned space"
    )]
    pub points_report: Option<PathBuf>,

    /// Contig to genome map in CAMI binning format, offered to the pool as extra candidates.
    /// A probe: it asks whether the bar would take the right grouping if handed one
    #[arg(long = "dissolve-oracle", hide_short_help = true)]
    pub dissolve_oracle: Option<String>,
}

#[derive(Args, Debug, Clone, Default)]
#[command(next_help_heading = "Reports")]
pub struct MarkerReport {
    /// Write every single copy marker hit, with whether its gene ran off a contig end
    #[arg(id = "marker_report", long = "marker-report", hide_short_help = true)]
    pub path: Option<PathBuf>,
}

// Probes of the first partition. A run given its bins has no first partition to probe.
#[derive(Args, Debug, Clone, Default)]
#[command(next_help_heading = "Reports")]
pub struct PartitionProbes {
    /// Write pieces of long contigs and real short contigs, per length band, against their
    /// nearest long contigs. Stops before the annotation
    #[arg(long = "floor-report", hide_short_help = true)]
    pub floor_report: Option<PathBuf>,

    /// Write every contig's nearest neighbours. Stops before refining
    #[arg(long = "knn-report", hide_short_help = true)]
    pub knn_report: Option<PathBuf>,

    /// Read a contig to genome map in CAMI binning format and write, per contig, how far the
    /// nearest contig of its own genome sits in the neighbour order. Stops before refining
    #[arg(long = "reach-report", num_args = 2, value_names = ["GOLD", "OUT"],
          hide_short_help = true)]
    pub reach_report: Option<Vec<PathBuf>>,

    /// Write every contig's community in every rung of every arm and seed, each arm's chosen
    /// rung and the combined bins, at both coverage weights. `PATH.pull` holds each contig's
    /// edge weight toward every combined community it touches. Stops before the outlier pool
    #[arg(long = "partition-report", hide_short_help = true)]
    pub partition_report: Option<PathBuf>,

    /// Write every candidate the ensemble offered the arbiter, with its worth and members.
    /// A probe: it asks whether a grouping the arbiter refused was ever proposed
    #[arg(long = "combine-report", hide_short_help = true)]
    pub combine_report: Option<PathBuf>,
}
