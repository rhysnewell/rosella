use std::env;

use clap::{Parser, crate_name, crate_version};
use env_logger::Builder;
use log::{LevelFilter, error, info};

use rosella::cli::{Cli, Command, Logging};
use rosella::pool;
use rosella::recover::recover_engine::run_recover;
use rosella::quality::bins::run_score;
use rosella::refine::refinery::run_refine;

fn main() {
    rosella::timing::start();

    match Cli::parse().command {
        Command::Recover(args) => {
            set_log_level(&args.logging);
            start_pool(args.common.threads);
            exit_on_error("Recover", pool::install(|| run_recover(*args)));
        }
        Command::Refine(args) => {
            set_log_level(&args.logging);
            start_pool(args.common.threads);
            exit_on_error("Refine", pool::install(|| run_refine(*args)));
        }
        Command::Score(args) => {
            set_log_level(&args.logging);
            start_pool(args.threads);
            exit_on_error("Score", pool::install(|| run_score(*args)));
        }
    }
}

fn start_pool(threads: usize) {
    if let Err(e) = pool::init(threads) {
        error!("Failed to build a thread pool of {threads}: {e}");
        std::process::exit(1);
    }
}

fn exit_on_error(subcommand: &str, outcome: anyhow::Result<()>) {
    if let Err(e) = outcome {
        error!("{} Failed with error: {}", subcommand, e);
        std::process::exit(1);
    }
}

fn set_log_level(logging: &Logging) {
    let log_level = if logging.quiet {
        LevelFilter::Error
    } else if logging.verbose {
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    };

    let mut builder = Builder::new();
    builder.filter_level(log_level);
    if let Ok(filters) = env::var("RUST_LOG") {
        builder.parse_filters(&filters);
    }
    if builder.try_init().is_err() {
        panic!("Failed to set log level - has it been specified multiple times?")
    }
    info!("{} version {}", crate_name!(), crate_version!());
}
