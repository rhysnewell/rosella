use clap::*;
use clap_complete::*;

mod common;
mod recover;
mod refine;

pub mod help;

pub use help::{recover_full_help, refine_full_help};

pub fn build_cli() -> Command {
    Command::new("rosella")
        .version(crate_version!())
        .author(crate::AUTHOR_AND_EMAIL)
        .about(format!(
            "Recover MAGs from contigs using UMAP and HDBSCAN clustering. (version {})",
            crate_version!()
        ))
        .arg_required_else_help(true)
        .subcommand(recover::command())
        .subcommand(refine::command())
        .subcommand(shell_completion())
}

fn shell_completion() -> Command {
    Command::new("shell-completion")
        .about("Generate a shell completion script for rosella")
        .arg(
            Arg::new("output-file")
                .short('o')
                .long("output-file")
                .required(true),
        )
        .arg(
            Arg::new("shell")
                .long("shell")
                .required(true)
                .value_parser(value_parser!(Shell)),
        )
        .args(common::logging_args())
}
