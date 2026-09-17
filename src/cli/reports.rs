use clap::Args;

/// None of these changes a bin. They write what a run decided, or hand it something to judge.
#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Reports")]
pub struct ReportPaths {
    /// Write every single copy marker hit, with whether its gene ran off a contig end
    #[arg(long = "marker-report", hide_short_help = true)]
    pub marker_report: Option<String>,

    /// Write every candidate the rescue pool judged, with its rank, verdict and members
    #[arg(long = "pool-report", hide_short_help = true)]
    pub pool_report: Option<String>,

    /// Write every contig's nearest neighbours to this path and stop before partitioning
    #[arg(long = "knn-report", hide_short_help = true)]
    pub knn_report: Option<std::path::PathBuf>,

    /// Write every contig against its neighbourhood, with its own bin's share, each rival's
    /// share, and what the recruitment claim makes of the pair. Ungated so the length can be swept
    #[arg(long = "audit-report", hide_short_help = true)]
    pub audit_report: Option<std::path::PathBuf>,

    /// Contig to genome map in CAMI binning format, offered to the pool as extra candidates.
    /// A probe: it asks whether the bar would take the right grouping if it were handed one
    #[arg(long = "dissolve-oracle", hide_short_help = true)]
    pub dissolve_oracle: Option<String>,
}
