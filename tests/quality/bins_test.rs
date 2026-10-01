use clap::Parser;
use rosella::cli::{Cli, Command};
use rosella::quality::bins::run_score;

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
