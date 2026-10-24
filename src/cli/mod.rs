use clap::parser::ValueSource;
use clap::{ArgMatches, Parser, Subcommand};

pub mod binning;
pub mod coverage;
pub mod ensemble;
pub mod make;
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
pub use ensemble::EnsembleArgs;
pub use make::MakeArgs;
pub use recover::{RecoverArgs, RunArgs};
pub use refine::RefineArgs;
pub use reports::{PartitionProbes, ReportPaths};
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
    /// Bin contigs into MAGs
    Recover(Box<RecoverArgs>),
    /// Refine bins from any binner through recover's refine cycle
    Refine(Box<RefineArgs>),
    /// Bin contigs with several binners' bins offered beside recover's own
    Ensemble(Box<EnsembleArgs>),
    /// Score bins on the single copy markers, with no gold standard
    Score(Box<ScoreArgs>),
    /// Write a depth directory from BAMs or reads, for recover and refine to read
    Make(Box<MakeArgs>),
}

pub fn kmer_size_ignored(matches: &ArgMatches) -> bool {
    matches.subcommand().is_some_and(|(_, given)| {
        given
            .try_contains_id("kmer_frequency_file")
            .unwrap_or(false)
            && given.value_source("kmer_size") == Some(ValueSource::CommandLine)
    })
}
