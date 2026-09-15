//! There were no CLI tests at all, and CI never invokes the binary. `debug_assert` is what
//! would have caught a flag declared without an action, and the required-argument matrix is
//! what the seven hand-repeated `required_unless_present_any` lists used to encode.

use clap::{CommandFactory, Parser};
use rosella::cli::Cli;
use rosella::clustering::graph_partition::{PARTITION_NAMES, Partition};

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

/// Supplying both tables does not excuse `refine` from `--assembly`: the bin writer reads
/// the sequence off it, so a run without one fails after the binning rather than before it.
#[test]
fn refine_needs_an_assembly_even_with_both_tables() {
    let error = parse(&[
        "refine", "-o", "out", "-f", "bin.fna", "-C", "cov.tsv", "-K", "kmers.tsv",
    ])
    .unwrap_err()
    .to_string();
    assert!(error.contains("--assembly"), "{error}");
}

#[test]
fn the_bounded_parsers_reject_values_outside_their_range() {
    for (flag, value) in [
        ("--knn-candidates", "1"),
        ("--knn-candidates", "257"),
        ("--split-level-quantile", "1.5"),
        ("--presence-fraction", "1.5"),
        ("--marker-bar-offset", "101"),
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

/// The engine expects clap to have restricted this and panics otherwise, so a name added to
/// one list and not the other is a crash rather than a rejected argument.
#[test]
fn every_partition_name_clap_accepts_has_a_parser_behind_it() {
    for name in PARTITION_NAMES {
        assert!(
            recover_with(&["-C", "cov.tsv", "--partition", name]).is_ok(),
            "clap rejected --partition {name}"
        );
        assert!(Partition::parse(name).is_some(), "nothing parses {name}");
    }
    assert!(recover_with(&["-C", "cov.tsv", "--partition", "nonesuch"]).is_err());
}


/// The splitter builds its per-bin graphs through the same features as `recover`, so the graph
/// has to reach both subcommands from one place.
#[test]
fn both_subcommands_take_an_assembly_graph() {
    let graph = ["--assembly-graph", "graph.gfa", "--assembly-graph-weight", "0.25"];
    for command in [
        vec!["recover", "-r", "a.fna", "-o", "out", "-C", "cov.tsv"],
        vec!["refine", "-r", "a.fna", "-o", "out", "-C", "cov.tsv", "-f", "bin.fna"],
    ] {
        let mut arguments = command.clone();
        arguments.extend_from_slice(&graph);
        let parsed = parse(&arguments).unwrap_or_else(|error| panic!("{command:?}: {error}"));
        let params = match parsed.command {
            rosella::cli::Command::Recover(args) => args.graph,
            rosella::cli::Command::Refine(args) => args.graph,
            _ => unreachable!(),
        };
        assert_eq!(params.assembly_graph.as_deref(), Some("graph.gfa"));
        assert_eq!(params.assembly_graph_weight, 0.25);
    }

    let mut weight_alone = vec!["refine", "-o", "out", "-f", "b.fna", "-C", "c.tsv", "-K", "k.tsv"];
    weight_alone.extend_from_slice(&graph[2..]);
    assert!(parse(&weight_alone).is_err(), "the weight needs a graph");
}
