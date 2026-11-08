use std::path::PathBuf;

use clap::Args;

#[derive(Args, Debug, Clone, Default)]
#[command(next_help_heading = "Reports")]
pub struct ReportPaths {
    #[command(flatten)]
    pub marker: MarkerReport,

    #[arg(long = "stage-report", hide = true)]
    pub stages: Option<PathBuf>,
}

#[derive(Args, Debug, Clone, Default)]
#[command(next_help_heading = "Reports")]
pub struct MarkerReport {
    /// Write every single copy marker hit, with whether its gene ran off a contig end
    #[arg(id = "marker_report", long = "marker-report", hide_short_help = true)]
    pub path: Option<PathBuf>,
}
