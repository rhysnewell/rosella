use clap::{Parser, Subcommand};
use clap_complete::Shell;

pub mod common;
pub mod manual;
pub mod recover;
pub mod refine;

pub use common::*;
pub use recover::RecoverArgs;
pub use refine::RefineArgs;

#[derive(Parser, Debug)]
#[command(
    name = "rosella",
    version = concat!(env!("CARGO_PKG_VERSION"), " (", env!("ROSELLA_BUILD_COMMIT"), ")"),
    author = crate::AUTHOR_AND_EMAIL,
    about = "Recover MAGs from contigs using composition and coverage.",
    arg_required_else_help = true
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
    /// Generate a shell completion script for rosella
    ShellCompletion(ShellCompletionArgs),
}

#[derive(clap::Args, Debug, Clone)]
pub struct ShellCompletionArgs {
    /// Where the completion script is written
    #[arg(short, long = "output-file")]
    pub output_file: String,

    /// Shell to generate for
    #[arg(long)]
    pub shell: Shell,

    #[command(flatten)]
    pub logging: Logging,
}
