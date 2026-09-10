use rosella::quality::cache::{Annotation, path_for, read, select, write};

fn names() -> Vec<String> {
    vec!["one".to_string(), "two".to_string()]
}

fn annotation() -> Annotation {
    Annotation {
        metadata: vec![[0u32; 22], {
            let mut row = [1u32; 22];
            row[20] = 4096;
            row[21] = 7;
            row
        }],
        hits: vec![Vec::new(), vec![(3, 1), (11, 2)]],
    }
}

#[test]
fn the_tables_come_back_as_they_went_in() {
    let home = tempfile::tempdir().unwrap();
    let path = home.path().join("families.gz");
    let held = annotation();
    write(&path, &names(), &held).unwrap();

    let (back, table) = read(&path).unwrap();
    assert_eq!(back, names());
    assert_eq!(table.metadata, held.metadata);
    assert_eq!(table.hits, held.hits);
}

#[test]
fn a_subset_comes_back_in_the_order_it_was_asked_for() {
    let held = annotation();
    let taken = select(&["two".to_string()], &names(), annotation()).unwrap();

    assert_eq!(taken.metadata, vec![held.metadata[1]]);
    assert_eq!(taken.hits, vec![held.hits[1].clone()]);
}

#[test]
fn a_contig_the_table_never_held_is_refused() {
    let wanted = ["one".to_string(), "three".to_string()];

    assert!(select(&wanted, &names(), annotation()).is_err());
}

#[test]
fn the_key_follows_the_assembly_the_database_and_the_tier() {
    let home = std::path::Path::new("/cache");
    let database = std::path::Path::new("uniref100.dmnd");

    let held = path_for(home, "assembly.fa", database, "default");
    assert_eq!(held, path_for(home, "assembly.fa", database, "default"));
    assert_ne!(held, path_for(home, "other.fa", database, "default"));
    assert_ne!(held, path_for(home, "assembly.fa", database, "faster"));
    assert_ne!(
        held,
        path_for(
            home,
            "assembly.fa",
            std::path::Path::new("other.dmnd"),
            "default"
        )
    );
}
