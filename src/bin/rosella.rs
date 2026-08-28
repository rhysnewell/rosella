use std::env;

use clap::{CommandFactory, Parser, crate_name, crate_version};
use clap_complete::generate;
use env_logger::Builder;
use log::{LevelFilter, error, info};

use rosella::cli::{Cli, Command, Logging, manual};
use rosella::recover::recover_engine::run_recover;
use rosella::refine::refinery::run_refine;

fn main() {
    rosella::timing::start();

    // Before clap parses, so a missing required argument cannot stop the manual printing.
    if let Some(request) = manual::requested() {
        if let Err(e) = manual::print(&request) {
            eprintln!("{e}");
            std::process::exit(1);
        }
        return;
    }

    match Cli::parse().command {
        Command::Recover(args) => {
            set_log_level(&args.logging);
            use_threads(args.common.threads);
            exit_on_error("Recover", run_recover(*args));
        }
        Command::Refine(args) => {
            set_log_level(&args.logging);
            use_threads(args.common.threads);
            exit_on_error("Refine", run_refine(*args));
        }
        Command::ShellCompletion(args) => {
            set_log_level(&args.logging);
            let mut file =
                std::fs::File::create(&args.output_file).expect("failed to open output file");
            let mut command = Cli::command();
            info!("Generating completion script for shell {}", args.shell);
            let name = command.get_name().to_string();
            generate(args.shell, &mut command, name, &mut file);
        }
    }
}

fn use_threads(threads: usize) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
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
