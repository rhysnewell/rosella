use std::io::Write;

use clap::Parser;
use rosella::cli::{Common, DistanceParams, ReadArgs};
use rosella::kmers::kmer_counting::prefixes;
use rosella::tables::{Lengths, Sources, Tables};

#[derive(Parser)]
struct Cli {
    #[command(flatten)]
    common: Common,
    #[command(flatten)]
    reads: ReadArgs,
    #[command(flatten)]
    distance: DistanceParams,
}

/// A contig of Ns leaves both tables by name and by row, so every row left still pairs a
/// contig's coverage with its own composition. A soft-masked contig is not blank.
#[test]
fn a_blank_contig_leaves_both_tables_aligned() {
    let home = tempfile::tempdir().unwrap();
    let assembly = home.path().join("assembly.fasta");
    let mut fasta = std::fs::File::create(&assembly).unwrap();
    for (name, sequence) in [
        ("masked", "acgttgcaaggcttac".repeat(100)),
        ("blank", "N".repeat(1600)),
        ("plain", "GGGCATTACCAGTTTA".repeat(100)),
    ] {
        writeln!(fasta, ">{name}\n{sequence}").unwrap();
    }
    let coverage = home.path().join("coverage.tsv");
    std::fs::write(
        &coverage,
        "contigName\tcontigLen\ttotalAvgDepth\ts1\ts1-var\n\
         masked\t1600\t1.0\t1.0\t0.1\n\
         blank\t1600\t2.0\t2.0\t0.1\n\
         plain\t1600\t3.0\t3.0\t0.1\n",
    )
    .unwrap();
    let output = home.path().join("out");
    let cli = Cli::parse_from([
        "rosella",
        "--output-directory",
        output.to_str().unwrap(),
        "--coverage-file",
        coverage.to_str().unwrap(),
    ]);

    let tables = Tables::build(&Sources {
        assembly: assembly.to_str().unwrap(),
        common: &cli.common,
        reads: &cli.reads,
        distance: &cli.distance,
        lengths: Lengths::from(0, 0),
        threads: 1,
    })
    .unwrap();

    assert_eq!(tables.coverage.contig_names, ["masked", "plain"]);
    assert_eq!(tables.coverage.average_depths, [1.0, 3.0]);
    let own = prefixes(
        assembly.to_str().unwrap(),
        &[("masked", 1600), ("plain", 1600)],
        cli.distance.kmer_size,
    )
    .unwrap();
    assert_eq!(tables.tnf.kmer_table, own);
}
