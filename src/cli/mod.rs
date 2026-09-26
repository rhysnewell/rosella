use clap::{Parser, Subcommand};

pub mod binning;
pub mod coverage;
pub mod markers;
pub mod recover;
pub mod refine;
pub mod reports;
pub mod rescue;
pub mod runtime;
pub mod score;
pub mod style;

pub use binning::*;
pub use coverage::*;
pub use markers::MarkerParams;
pub use recover::RecoverArgs;
pub use refine::RefineArgs;
pub use reports::ReportPaths;
pub use rescue::RescueParams;
pub use runtime::*;
pub use score::ScoreArgs;

#[derive(Parser, Debug)]
#[command(
    name = "rosella",
    version = concat!(env!("CARGO_PKG_VERSION"), " (", env!("ROSELLA_BUILD_COMMIT"), ")"),
    author = crate::AUTHOR_AND_EMAIL,
    about = "Recover MAGs from contigs using composition and coverage.",
    arg_required_else_help = true,
    styles = style::HELP
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Recover MAGs from contigs by partitioning the k-nearest-neighbour graph.
    Recover(Box<RecoverArgs>),
    /// Refine MAGs by re-partitioning each bin on its own.
    Refine(Box<RefineArgs>),
    /// Score a set of bins against the single copy markers, with no gold standard
    Score(Box<ScoreArgs>),
}
