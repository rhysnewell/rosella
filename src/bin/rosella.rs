use clap::crate_version;
use env_logger::Builder;
use log::{LevelFilter, info, error};
use std::env;

use rosella::cli::build_cli;
use rosella::recover::recover_engine::run_recover;

fn main() {
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    set_log_level(&matches, false);

    match matches.subcommand_name() {
        Some("recover") => {
            let sub_matches = matches.subcommand_matches("recover").unwrap();
            set_log_level(&sub_matches, true);
            match run_recover(sub_matches) {
                Ok(_) => {}
                Err(e) => {
                    error!("Recover Failed with error: {}", e);
                    std::process::exit(1);
                }
            };
        },
        Some("refine") => {
            unimplemented!();
        },
        _ => {
            app.print_help().unwrap();
            std::process::exit(1);
        }
    }
}


fn set_log_level(matches: &clap::ArgMatches, is_last: bool) {
    let mut log_level = LevelFilter::Info;
    let mut specified = false;
    if matches.get_flag("verbose") {
        specified = true;
        log_level = LevelFilter::Debug;
    }
    if matches.get_flag("quiet") {
        specified = true;
        log_level = LevelFilter::Error;
    }
    if specified || is_last {
        let mut builder = Builder::new();
        builder.filter_level(log_level);
        if env::var("RUST_LOG").is_ok() {
            builder.parse_filters(&env::var("RUST_LOG").unwrap());
        }
        if builder.try_init().is_err() {
            panic!("Failed to set log level - has it been specified multiple times?")
        }
    }
    if is_last {
        info!("lorikeet version {}", crate_version!());
    }
}