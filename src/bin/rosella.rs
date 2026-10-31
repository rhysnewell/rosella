use clap::{ArgMatches, CommandFactory, FromArgMatches, crate_name, crate_version};
use log::{LevelFilter, error, info, warn};

use rosella::cli::{Cli, Command, Logging, kmer_size_ignored};
use rosella::coverage::make::run_make;
use rosella::pool;
use rosella::quality::bins::run_score;
use rosella::recover::recover_engine::{run_ensemble, run_recover};
use rosella::refine::refinery::run_refine;

fn main() {
    rosella::timing::start();

    let matches = Cli::command().get_matches();
    let cli = Cli::from_arg_matches(&matches).unwrap_or_else(|error| error.exit());
    match cli.command {
        Command::Recover(args) => run(
            "Recover",
            &matches,
            &args.run.logging,
            args.run.runtime.threads,
            || run_recover(&args),
        ),
        Command::Refine(args) => run(
            "Refine",
            &matches,
            &args.run.logging,
            args.run.runtime.threads,
            || run_refine(&args),
        ),
        Command::Ensemble(args) => run(
            "Ensemble",
            &matches,
            &args.recover.run.logging,
            args.recover.run.runtime.threads,
            || run_ensemble(&args),
        ),
        Command::Score(args) => run(
            "Score",
            &matches,
            &args.logging,
            args.runtime.threads,
            || run_score(&args),
        ),
        Command::Make(args) => run(
            "Make",
            &matches,
            &args.logging,
            args.runtime.threads,
            || run_make(&args),
        ),
    }
}

fn run(
    subcommand: &str,
    matches: &ArgMatches,
    logging: &Logging,
    threads: usize,
    job: impl FnOnce() -> anyhow::Result<()> + Send,
) {
    set_log_level(logging);
    if kmer_size_ignored(matches) {
        warn!("--kmer-size: ignored beside -K");
    }
    if let Err(e) = pool::init(threads) {
        error!("Thread pool of {threads}: {e}");
        std::process::exit(1);
    }
    if let Err(e) = pool::install(job) {
        error!("{subcommand} failed: {e:#}");
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
    info!("Version: {} {}", crate_name!(), crate_version!());
}
