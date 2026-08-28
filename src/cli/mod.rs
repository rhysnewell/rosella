use clap::{Parser, Subcommand, crate_version};
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
    version = crate_version!(),
    author = crate::AUTHOR_AND_EMAIL,
    about = "Recover MAGs from contigs using UMAP and HDBSCAN clustering.",
    arg_required_else_help = true
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Recover MAGs from contigs using UMAP and HDBSCAN clustering.
    Recover(Box<RecoverArgs>),
    /// Refine MAGs using UMAP and HDBSCAN clustering.
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
