use clap::{crate_version, crate_name};
use env_logger::Builder;
use log::{LevelFilter, info, error};
use rosella::refine::refinery::run_refine;
use std::env;

use rosella::cli::{build_cli, refine_full_help};
use rosella::recover::recover_engine::run_recover;

use bird_tool_utils::clap_utils::print_full_help_if_needed;

fn main() {
    let mut app = build_cli();
    let matches = app.clone().get_matches();

    match matches.subcommand_name() {
        Some("recover") => {
            let sub_matches = matches.subcommand_matches("recover").unwrap();
            print_full_help_if_needed(sub_matches, refine_full_help());
            set_log_level(&sub_matches, true);
            // set rayon threads
            let threads = *sub_matches.get_one::<usize>("threads").unwrap();
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
            match run_recover(sub_matches) {
                Ok(_) => {}
                Err(e) => {
                    error!("Recover Failed with error: {}", e);
                    std::process::exit(1);
                }
            };
        },
        Some("refine") => {
            #[cfg(feature = "no_flight")]
            {
                let sub_matches = matches.subcommand_matches("refine").unwrap();
                print_full_help_if_needed(sub_matches, refine_full_help());
                set_log_level(&sub_matches, true);
                // set rayon threads
                let threads = *sub_matches.get_one::<usize>("threads").unwrap();
                rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
                error!("Refine is not available in this version of rosella");
                error!("Recompile without the 'no_flight' feature and install flight via GitHub");
                unimplemented!();
            }

            #[cfg(not(feature = "no_flight"))]
            {
                let sub_matches = matches.subcommand_matches("refine").unwrap();
                print_full_help_if_needed(sub_matches, refine_full_help());
                set_log_level(&sub_matches, true);
                // set rayon threads
                let threads = *sub_matches.get_one::<usize>("threads").unwrap();
                rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
                match run_refine(sub_matches) {
                    Ok(_) => {}
                    Err(e) => {
                        error!("Refine Failed with error: {}", e);
                        std::process::exit(1);
                    }
                };
            }
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
        builder.filter_module("annembed", LevelFilter::Off);
        builder.filter_module("hnsw_rs", LevelFilter::Off);
        if env::var("RUST_LOG").is_ok() {
            builder.parse_filters(&env::var("RUST_LOG").unwrap());
        }
        if builder.try_init().is_err() {
            panic!("Failed to set log level - has it been specified multiple times?")
        }
    }
    if is_last {
        info!("{} version {}", crate_name!(), crate_version!());
    }
}