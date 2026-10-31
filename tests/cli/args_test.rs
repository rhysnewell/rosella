//! There were no CLI tests at all, and CI never invokes the binary. `debug_assert` is what
//! would have caught a flag declared without an action, and the required-argument matrix is
//! what the seven hand-repeated `required_unless_present_any` lists used to encode.

use clap::{CommandFactory, Parser};
use rosella::cli::{Cli, kmer_size_ignored};

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
    .unwrap_err()
    .to_string();
    assert!(error.contains("--assembly"), "{error}");
}

#[test]
fn the_bounded_parsers_reject_values_outside_their_range() {
    for (flag, value) in [
        ("--min-completeness", "101"),
        ("--marker-fragment-span", "1.5"),
    ] {
        let error = recover_with(&["-C", "cov.tsv", flag, value])
            .unwrap_err()
            .to_string();
        assert!(error.contains("outside"), "{flag} {value}: {error}");
    }
}

/// The splitter builds its per-bin graphs through the same features as `recover`, so the graph
/// has to reach both subcommands from one place.
#[test]
fn both_subcommands_take_an_assembly_graph() {
    let graph = [
        "--assembly-graph",
        "graph.gfa",
        "--assembly-graph-weight",
        "0.25",
    ];
    for command in [
        vec!["recover", "-r", "a.fna", "-o", "out", "-C", "cov.tsv"],
        vec![
            "refine", "-r", "a.fna", "-o", "out", "-C", "cov.tsv", "-f", "bin.fna",
        ],
    ] {
        let mut arguments = command.clone();
        arguments.extend_from_slice(&graph);
        let parsed = parse(&arguments).unwrap_or_else(|error| panic!("{command:?}: {error}"));
        let params = match parsed.command {
            rosella::cli::Command::Recover(args) => args.run.graph,
            rosella::cli::Command::Refine(args) => args.run.graph,
            _ => unreachable!(),
        };
        assert_eq!(params.assembly_graph.as_deref(), Some("graph.gfa"));
        assert_eq!(params.assembly_graph_weight, 0.25);
    }

    let mut weight_alone = vec![
        "refine", "-o", "out", "-f", "b.fna", "-C", "c.tsv", "-K", "k.tsv",
    ];
    weight_alone.extend_from_slice(&graph[2..]);
    assert!(parse(&weight_alone).is_err(), "the weight needs a graph");
}

/// refine starts from bins, so it reads every run flag recover does but has no first partition
/// to seed or probe.
#[test]
fn refine_takes_the_run_flags_but_not_the_partition_probes() {
    let refine = [
        "refine", "-r", "a.fna", "-o", "out", "-C", "cov.tsv", "-f", "bin.fna",
    ];
    let shared = ["--marker-cache", "cache", "--no-refine", "--gpu"];
    let parsed = parse(&[refine.as_slice(), &shared].concat()).expect("refine parses");
    let rosella::cli::Command::Refine(args) = parsed.command else {
        unreachable!()
    };
    assert_eq!(args.run.markers.marker_cache.as_deref(), Some("cache"));
    assert!(args.run.rescue.no_refine && args.run.binning.gpu);
    for probe in [["--knn-report", "k.tsv"], ["--partition-seeds", "2"]] {
        assert!(
            parse(&[refine.as_slice(), &probe].concat()).is_err(),
            "{probe:?}"
        );
    }
}

/// make writes the depth directory, so a table or a directory to read from has no place in it.
#[test]
fn make_takes_reads_and_bams_but_no_coverage_table() {
    assert!(parse(&["make", "-o", "depth", "-b", "a.bam"]).is_ok());
    assert!(parse(&["make", "-o", "depth"]).is_ok());
    assert!(parse(&["make", "-o", "depth", "-C", "cov.tsv"]).is_err());
}

#[test]
fn an_attach_floor_never_reaches_past_the_min_contig_size() {
    for (flags, cutoff, floor) in [
        (vec![], 1500, None),
        (vec!["--attach-floor", "500"], 1500, Some(500)),
        (
            vec!["--min-contig-size", "1000", "--attach-floor", "2000"],
            1000,
            Some(1000),
        ),
        (vec!["--min-contig-size", "600"], 600, None),
    ] {
        let parsed = recover_with(&[["-C", "cov.tsv"].as_slice(), &flags].concat())
            .unwrap_or_else(|error| panic!("{flags:?}: {error}"));
        let rosella::cli::Command::Recover(args) = parsed.command else {
            unreachable!()
        };
        assert_eq!(
            (args.run.binning.cutoff(), args.run.binning.attach_floor()),
            (cutoff, floor),
            "{flags:?}"
        );
    }
    assert!(recover_with(&["-C", "cov.tsv", "--attach-floor", "short"]).is_err());
}

#[test]
fn a_kmer_size_typed_beside_a_kmer_table_is_flagged() {
    let matches = |extra: &[&str]| {
        let mut arguments = vec![
            "rosella", "recover", "-r", "a.fna", "-o", "out", "-C", "c.tsv",
        ];
        arguments.extend_from_slice(extra);
        Cli::command().try_get_matches_from(arguments).unwrap()
    };

    assert!(kmer_size_ignored(&matches(&[
        "-K",
        "k.tsv",
        "--kmer-size",
        "5"
    ])));
    assert!(!kmer_size_ignored(&matches(&["-K", "k.tsv"])));
    assert!(!kmer_size_ignored(&matches(&["--kmer-size", "5"])));
}
