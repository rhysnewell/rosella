use std::fs::File;

use rosella::bins::{discover, refuse_used};
use tempfile::tempdir;

fn touch(directory: &std::path::Path, name: &str) -> String {
    let path = directory.join(name);
    File::create(&path).unwrap();
    path.to_string_lossy().into_owned()
}

#[test]
fn a_leading_dot_on_the_extension_still_matches() {
    let home = tempdir().unwrap();
    touch(home.path(), "one.fna");
    touch(home.path(), "two.fna");
    touch(home.path(), "notes.txt");

    let directory = home.path().to_string_lossy().into_owned();
    let dotted = discover(&[], Some(&directory), ".fna").unwrap();
    let plain = discover(&[], Some(&directory), "fna").unwrap();
    assert_eq!(dotted, plain);
    assert_eq!(dotted.len(), 2);
}

#[test]
fn files_and_directory_are_read_together_and_deduplicated() {
    let home = tempdir().unwrap();
    let one = touch(home.path(), "one.fna");
    touch(home.path(), "two.fna");

    let directory = home.path().to_string_lossy().into_owned();
    let found = discover(&[one], Some(&directory), "fna").unwrap();
    assert_eq!(found.len(), 2, "{found:?}");
    assert!(found.windows(2).all(|pair| pair[0] < pair[1]));
}

#[test]
fn naming_nothing_is_an_error_rather_than_an_empty_run() {
    let home = tempdir().unwrap();
    let directory = home.path().to_string_lossy().into_owned();
    assert!(discover(&[], Some(&directory), "fna").is_err());
}

#[test]
fn a_directory_holding_bins_refuses_a_second_run() {
    let home = tempdir().unwrap();
    let directory = home.path().to_string_lossy().into_owned();
    assert!(refuse_used(&directory).is_ok());
    touch(home.path(), "rosella_bin_1.fna");
    assert!(refuse_used(&directory).is_err());
    assert!(refuse_used("a directory that does not exist").is_ok());
}
