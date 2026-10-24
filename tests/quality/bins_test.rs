use clap::Parser;
use rosella::cli::{Cli, Command};
use std::collections::HashMap;

use rosella::quality::bins::run_score;
use rosella::quality::orfs::call_over;

/// Bins are keyed by file stem, so two bins with one stem used to leave one row for both.
#[test]
fn two_bins_with_one_name_are_refused() {
    let home = tempfile::tempdir().unwrap();
    let mut paths = Vec::new();
    for (folder, contig) in [("a", "c1"), ("b", "c2")] {
        let path = home.path().join(folder).join("bin.fna");
        std::fs::create_dir(path.parent().unwrap()).unwrap();
        std::fs::write(&path, format!(">{contig}\nACGT\n")).unwrap();
        paths.push(path.to_string_lossy().into_owned());
    }
    let output = home.path().join("quality.tsv");
    let Command::Score(args) = Cli::parse_from([
        "rosella",
        "score",
        "-f",
        &paths[0],
        &paths[1],
        "-o",
        &output.to_string_lossy(),
    ])
    .command
    else {
        unreachable!("score parses as score");
    };

    let error = run_score(&args).unwrap_err().to_string();
    assert!(error.contains("both named bin"), "{error}");
}

#[test]
fn score_writes_every_called_protein_numbered_within_its_contig() {
    rosella::pool::init(2).ok();
    let home = tempfile::tempdir().unwrap();
    let assembly = home.path().join("joined.fna");
    let mut pieces = Vec::new();
    let mut reader = needletail::parse_fastx_file("tests/data/ben/1.fna.gz").unwrap();
    while let Some(record) = reader.next()
        && pieces.len() < 400
    {
        pieces.push(record.unwrap().seq().into_owned());
    }
    let (first, second) = pieces.split_at(200);
    std::fs::write(
        &assembly,
        [
            b">a\n".as_slice(),
            &first.concat(),
            b"\n>b\n",
            &second.concat(),
            b"\n",
        ]
        .concat(),
    )
    .unwrap();
    let output = home.path().join("quality.tsv");
    let proteins = home.path().join("proteins.faa");
    let Command::Score(args) = Cli::parse_from([
        "rosella",
        "score",
        "-f",
        &assembly.to_string_lossy(),
        "-o",
        &output.to_string_lossy(),
        "--proteins",
        &proteins.to_string_lossy(),
    ])
    .command
    else {
        unreachable!("score parses as score");
    };
    run_score(&args).unwrap();

    let mut called = Vec::new();
    call_over(
        &assembly.to_string_lossy(),
        0..usize::MAX,
        |batch, names| {
            called.extend(
                batch
                    .into_iter()
                    .map(|orf| (names[orf.contig].clone(), orf.protein)),
            );
            Ok(())
        },
    )
    .unwrap();
    let written = std::fs::read_to_string(&proteins).unwrap();
    let lines = written.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 2 * called.len());
    let mut genes = HashMap::new();
    for (record, (contig, protein)) in lines.chunks(2).zip(&called) {
        let gene = genes.entry(contig.as_str()).or_insert(0);
        *gene += 1;
        assert_eq!(
            record,
            [format!(">{contig}_{gene}").as_str(), protein.as_str()]
        );
    }
    assert!(genes.values().all(|count| *count > 1), "{genes:?}");
}
