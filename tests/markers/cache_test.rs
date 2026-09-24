use std::fs;
use std::os::unix::fs::symlink;
use std::path::Path;
use std::thread;

use rosella::markers::cache::{find, key, read, write};
use rosella::markers::replicon::Shape;
use rosella::markers::{Hit, MarkerSet};
fn entry(directory: &Path, header: &str) {
    fs::write(
        directory.join("markers.0123456789abcdef.tsv"),
        format!("rosella-markers-3\t{header}\ncontig_1\tPF00001:0\t9000\t10\n"),
    )
    .unwrap();
}

fn spelled(key: &str, path: &str) -> String {
    let mut fields = key.split('\t').collect::<Vec<_>>();
    fields[3] = path;
    fields.join("\t")
}

#[test]
fn an_entry_written_under_another_spelling_is_still_found() {
    let home = tempfile::tempdir().unwrap();
    let assembly = home.path().join("assembly.fasta");
    fs::write(&assembly, ">contig_1\nACGT\n").unwrap();
    let link = home.path().join("link.fasta");
    symlink(&assembly, &link).unwrap();

    let cache = home.path().join("cache");
    fs::create_dir(&cache).unwrap();

    let wanted = key(assembly.to_str().unwrap(), 1500, 0.3).unwrap();
    entry(&cache, &spelled(&wanted, link.to_str().unwrap()));

    assert!(find(&cache, &wanted).is_some());
    let through_link = key(link.to_str().unwrap(), 1500, 0.3).unwrap();
    assert!(find(&cache, &through_link).is_some());
}

#[test]
fn an_entry_taken_under_other_settings_is_refused() {
    let home = tempfile::tempdir().unwrap();
    let assembly = home.path().join("assembly.fasta");
    fs::write(&assembly, ">contig_1\nACGT\n").unwrap();

    let cache = home.path().join("cache");
    fs::create_dir(&cache).unwrap();

    let assembly = assembly.to_str().unwrap();
    entry(&cache, &key(assembly, 1500, 0.3).unwrap());

    assert!(find(&cache, &key(assembly, 1500, 0.5).unwrap()).is_none());
    assert!(find(&cache, &key(assembly, 2500, 0.3).unwrap()).is_none());
}

#[test]
fn a_directory_of_other_files_holds_nothing() {
    let home = tempfile::tempdir().unwrap();
    let assembly = home.path().join("assembly.fasta");
    fs::write(&assembly, ">contig_1\nACGT\n").unwrap();

    let cache = home.path().join("markers");
    fs::create_dir(&cache).unwrap();
    fs::write(cache.join("bacteria.hmm"), "HMMER3/f\n").unwrap();

    let wanted = key(assembly.to_str().unwrap(), 1500, 0.3).unwrap();
    assert!(find(&cache, &wanted).is_none());
}

#[test]
fn two_writers_on_one_entry_leave_one_whole_file() {
    let home = tempfile::tempdir().unwrap();
    let path = home.path().join("markers.0123456789abcdef.tsv");
    let table = "model_name\tdomain\nalpha\tbac120\n";
    let contigs = 20_000;

    thread::scope(|scope| {
        for writer in ["a", "b"] {
            let path = &path;
            scope.spawn(move || {
                let set = MarkerSet::parse(table);
                let names = (0..contigs).map(|at| format!("{writer}{at}")).collect::<Vec<_>>();
                let hits = vec![vec![Hit { marker: 0, partial: false }]; contigs];
                let shapes = vec![Shape { coding_bases: 900, genes: 1 }; contigs];
                for _ in 0..5 {
                    write(path, "key", &set, &names, &hits, &shapes).unwrap();
                }
            });
        }
    });

    let (names, _, _) = read(&path, &MarkerSet::parse(table)).unwrap();
    assert_eq!(names.len(), contigs);
    let writer = &names[0][..1];
    assert!(names.iter().all(|name| name.starts_with(writer)));
    assert_eq!(fs::read_dir(home.path()).unwrap().count(), 1);
}
