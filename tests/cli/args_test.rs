//! There were no CLI tests at all, and CI never invokes the binary. `debug_assert` is what
//! would have caught a flag declared without an action, and the required-argument matrix is
//! what the seven hand-repeated `required_unless_present_any` lists used to encode.

use clap::{CommandFactory, Parser};
use rosella::cli::{Cli, Command, manual};
use rosella::clustering::graph_partition::{NODE_SIZE_NAMES, NodeSize, PARTITION_NAMES, Partition};
use rosella::clustering::objective::{OBJECTIVE_NAMES, ObjectiveChoice};
use rosella::recover::ladder::{COMBINE_SOURCE_NAMES, CombineSource, RUNG_STATISTIC_NAMES, RungStatistic};

fn parse(arguments: &[&str]) -> Result<Cli, clap::Error> {
    Cli::try_parse_from(std::iter::once("rosella").chain(arguments.iter().copied()))
}

fn recover_with(source: &[&str]) -> Result<Cli, clap::Error> {
    let mut arguments = vec!["recover", "-r", "assembly.fna", "-o", "out"];
    arguments.extend_from_slice(source);
    parse(&arguments)
}

#[test]
fn the_command_tree_is_well_formed() {
    Cli::command().debug_assert();
}

#[test]
fn any_single_coverage_source_is_enough() {
    for source in [
        vec!["--coupled", "a.fq", "b.fq"],
        vec!["--interleaved", "a.fq"],
        vec!["--single", "a.fq"],
        vec!["--longreads", "a.fq"],
        vec!["--bam-files", "a.bam"],
        vec!["--longread-bam-files", "a.bam"],
        vec!["--coverage-file", "coverage.tsv"],
        vec!["--read1", "a.fq", "--read2", "b.fq"],
    ] {
        assert!(recover_with(&source).is_ok(), "{source:?} was rejected");
    }
}

#[test]
fn a_run_with_no_coverage_source_is_rejected() {
    let error = recover_with(&[]).unwrap_err().to_string();
    assert!(error.contains("--coverage-file"), "{error}");
}

/// The two satisfy each other rather than the group, so one alone is not a coverage source.
#[test]
fn forward_reads_alone_are_rejected() {
    let error = recover_with(&["--read1", "a.fq"]).unwrap_err().to_string();
    assert!(error.contains("--read2"), "{error}");
}

#[test]
fn refine_needs_bins_to_refine() {
    let error = parse(&["refine", "-r", "a.fna", "-o", "out", "-C", "cov.tsv"])
        .unwrap_err()
        .to_string();
    assert!(error.contains("--genome-fasta-files"), "{error}");
}

/// `refine` can work from the two tables alone, without ever reading the assembly.
#[test]
fn refine_takes_both_tables_in_place_of_an_assembly() {
    assert!(
        parse(&[
            "refine",
            "-o",
            "out",
            "-f",
            "bin.fna",
            "-C",
            "cov.tsv",
            "-K",
            "kmers.tsv",
        ])
        .is_ok()
    );
    let error = parse(&["refine", "-o", "out", "-f", "bin.fna", "-C", "cov.tsv"])
        .unwrap_err()
        .to_string();
    assert!(error.contains("--assembly"), "{error}");
}

#[test]
fn the_embedding_overrides_reject_values_outside_their_range() {
    for (flag, value) in [
        ("--umap-b", "9.0"),
        ("--umap-a", "0.0"),
        ("--min-dist", "6.0"),
        ("--spread", "0.0"),
    ] {
        let error = recover_with(&["-C", "cov.tsv", flag, value])
            .unwrap_err()
            .to_string();
        assert!(error.contains("outside"), "{flag} {value}: {error}");
    }
}

/// The mapper name is CoverM's to validate, so rosella must not reject one it has not
/// heard of.
#[test]
fn any_mapper_name_is_accepted() {
    for mapper in ["strobealign", "rammap-ont", "something-coverm-added-later"] {
        assert!(recover_with(&["-C", "cov.tsv", "-p", mapper]).is_ok());
    }
}

/// The engine expects clap to have restricted these and panics otherwise, so a name added to
/// one list and not the other is a crash rather than a rejected argument.
#[test]
fn every_name_clap_accepts_has_a_parser_behind_it() {
    for (flag, names, parse_one) in [
        (
            "--partition",
            PARTITION_NAMES.as_slice(),
            &(|name: &str| Partition::parse(name).is_some()) as &dyn Fn(&str) -> bool,
        ),
        (
            "--node-size",
            NODE_SIZE_NAMES.as_slice(),
            &(|name: &str| NodeSize::parse(name).is_some()),
        ),
        (
            "--objective",
            OBJECTIVE_NAMES.as_slice(),
            &(|name: &str| ObjectiveChoice::parse(name).is_some()),
        ),
        (
            "--rung-statistic",
            RUNG_STATISTIC_NAMES.as_slice(),
            &(|name: &str| RungStatistic::parse(name).is_some()),
        ),
    ] {
        for name in names {
            assert!(
                recover_with(&["-C", "cov.tsv", flag, name]).is_ok(),
                "clap rejected {flag} {name}"
            );
            assert!(parse_one(name), "nothing parses {flag} {name}");
        }
        assert!(recover_with(&["-C", "cov.tsv", flag, "nonesuch"]).is_err(), "{flag}");
    }
}

#[test]
fn both_subcommands_render_a_manual() {
    for subcommand in ["recover", "refine"] {
        let roff = manual::render(subcommand).unwrap();
        let roff = String::from_utf8(roff).unwrap();
        assert!(
            roff.contains(&format!("rosella-{subcommand}")),
            "{subcommand}"
        );
        assert!(roff.contains(".SH OPTIONS"), "{subcommand}");
    }
    assert!(manual::render("nonesuch").is_err());
}

#[test]
fn shell_completion_still_parses() {
    let cli = parse(&["shell-completion", "-o", "out.bash", "--shell", "bash"]).unwrap();
    assert!(matches!(cli.command, Command::ShellCompletion(_)));
}

/// Combining is the markers judging the ladder per bin, so turning the judge off silently took
/// `--combine-bins` with it.
#[test]
fn combining_bins_cannot_be_asked_for_with_the_marker_judge_off() {
    assert!(recover_with(&["-C", "cov.tsv", "--combine-bins"]).is_ok());
    assert!(recover_with(&["-C", "cov.tsv", "--no-marker-rungs"]).is_ok());
    assert!(recover_with(&["-C", "cov.tsv", "--no-marker-rungs", "--combine-bins"]).is_err());
}

/// `--combine-source` only means anything with combining on, so it carries `requires` and cannot
/// join the loop above, but it still needs every name it accepts to have a parser.
#[test]
fn every_combine_source_clap_accepts_has_a_parser_behind_it() {
    for name in COMBINE_SOURCE_NAMES {
        assert!(
            recover_with(&["-C", "cov.tsv", "--combine-bins", "--combine-source", name]).is_ok(),
            "clap rejected --combine-source {name}"
        );
        assert!(CombineSource::parse(name).is_some(), "nothing parses {name}");
    }
    assert!(recover_with(&["-C", "cov.tsv", "--combine-source", "arms"]).is_err());
}
