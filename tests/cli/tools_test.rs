use std::process::Command;

#[test]
fn a_missing_hmmsearch_stops_recover_before_it_reads_anything() {
    let scratch = tempfile::tempdir().unwrap();
    let output = scratch.path().join("out");
    let run = Command::new(env!("CARGO_BIN_EXE_rosella"))
        .args(["recover", "-r", "tests/data/ben/1.fna.gz"])
        .args(["--coverage-file", "absent.tsv", "-o"])
        .arg(&output)
        .env("PATH", "")
        .output()
        .unwrap();
    let said = String::from_utf8_lossy(&run.stderr);
    assert!(!run.status.success());
    assert!(said.contains("hmmsearch is not on PATH"), "{said}");
    assert!(!output.exists());
}
