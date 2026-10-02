use clap::{ArgMatches, CommandFactory, FromArgMatches, crate_name, crate_version};
use log::{LevelFilter, error, info, warn};

use rosella::cli::{Cli, Command, Logging, kmer_size_ignored};
use rosella::pool;
use rosella::quality::bins::run_score;
use rosella::recover::recover_engine::run_recover;
use rosella::refine::refinery::run_refine;

fn main() {
    rosella::timing::start();

    let matches = Cli::command().get_matches();
    let cli = Cli::from_arg_matches(&matches).unwrap_or_else(|error| error.exit());
    match cli.command {
        Command::Recover(args) => {
            set_log_level(&args.logging);
            warn_unread(&matches);
            start_pool(args.runtime.threads);
            exit_on_error("Recover", pool::install(|| run_recover(&args)));
        }
        Command::Refine(args) => {
            set_log_level(&args.logging);
            warn_unread(&matches);
            start_pool(args.runtime.threads);
            exit_on_error("Refine", pool::install(|| run_refine(&args)));
        }
        Command::Score(args) => {
            set_log_level(&args.logging);
            start_pool(args.runtime.threads);
            exit_on_error("Score", pool::install(|| run_score(&args)));
        }
    }
}

fn warn_unread(matches: &ArgMatches) {
    if kmer_size_ignored(matches) {
        warn!("--kmer-size is not read beside -K, whose table sets the k-mer size");
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
        error!("{subcommand} Failed with error: {e:#}");
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

    if rosella::progress::install(log_level).is_err() {
        panic!("Failed to set log level - has it been specified multiple times?")
    }
    info!("{} version {}", crate_name!(), crate_version!());
}
