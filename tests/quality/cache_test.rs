use rosella::quality::cache::{Annotation, path_for, read, write};

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

    let back = read(&path, &names()).unwrap();
    assert_eq!(back.metadata, held.metadata);
    assert_eq!(back.hits, held.hits);
}

#[test]
fn a_subset_comes_back_in_the_order_it_was_asked_for() {
    let home = tempfile::tempdir().unwrap();
    let path = home.path().join("families.gz");
    let held = annotation();
    write(&path, &names(), &held).unwrap();

    let back = read(&path, &["two".to_string()]).unwrap();
    assert_eq!(back.metadata, vec![held.metadata[1]]);
    assert_eq!(back.hits, vec![held.hits[1].clone()]);
}

#[test]
fn a_contig_the_table_never_held_is_refused() {
    let home = tempfile::tempdir().unwrap();
    let path = home.path().join("families.gz");
    write(&path, &names(), &annotation()).unwrap();

    assert!(read(&path, &["one".to_string(), "three".to_string()]).is_err());
}

#[test]
fn the_key_follows_the_assembly_and_the_database() {
    let home = std::path::Path::new("/cache");
    let database = std::path::Path::new("uniref100.dmnd");

    let held = path_for(home, "assembly.fa", database);
    assert_eq!(held, path_for(home, "assembly.fa", database));
    assert_ne!(held, path_for(home, "other.fa", database));
    assert_ne!(
        held,
        path_for(home, "assembly.fa", std::path::Path::new("other.dmnd"))
    );
}
