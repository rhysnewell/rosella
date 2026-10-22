use std::collections::{HashMap, HashSet};
use std::io::Write;

use rosella::recover::origin::Given;

/// Bins are read by contig name and placed on the rows the tables hold. A contig the tables left
/// out drops from its bin, a contig in two bins stays in the first, and the quality table follows
/// its bin by file stem.
#[test]
fn given_bins_land_on_the_rows_the_tables_hold() {
    let home = tempfile::tempdir().unwrap();
    let bin = |name: &str, contigs: &[&str]| {
        let path = home.path().join(name);
        let mut file = std::fs::File::create(&path).unwrap();
        for contig in contigs {
            writeln!(file, ">{contig} described\nACGT").unwrap();
        }
        path
    };
    let paths = [
        bin("a.fna", &["one", "two", "gone"]),
        bin("b.fna", &["two", "three"]),
    ];
    let quality = HashMap::from([("b".to_string(), 20.0)]);
    let given = Given::read(&paths, "rosella_t_".into(), quality, Some(15.0)).unwrap();

    let placed = given.place(&["three", "one", "two"].map(String::from), HashSet::new());
    assert_eq!(placed.bins.get(&0), Some(&vec![1, 2]));
    assert_eq!(placed.bins.get(&1), Some(&vec![0]));
    assert_eq!(placed.contamination, HashMap::from([(1, 20.0)]));
}

#[test]
fn a_cami_binning_reads_as_bins_in_first_seen_order() {
    let home = tempfile::tempdir().unwrap();
    let path = home.path().join("binner.binning");
    std::fs::write(
        &path,
        "@Version:0.9.1\n@SampleID:s\n\n@@SEQUENCEID\tBINID\n\
         two\tz\none described\ta\nthree\tz\ntwo\ta\ngone\ta\n",
    )
    .unwrap();
    let given = Given::table(&path).unwrap();

    assert_eq!(given.len(), 2);
    let placed = given.place(&["three", "one", "two"].map(String::from), HashSet::new());
    assert_eq!(placed.bins.get(&0), Some(&vec![0, 2]));
    assert_eq!(placed.bins.get(&1), Some(&vec![1]));
}
