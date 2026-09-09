//! The gene family tables, written once and read back instead of searched again.

use rosella::quality::cache::{Annotation, path_for, read, write};

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
    write(&path, &held).unwrap();

    let back = read(&path, 2).unwrap();
    assert_eq!(back.metadata, held.metadata);
    assert_eq!(back.hits, held.hits);
}

/// Reusing a table written for a different contig set would index gene families by the wrong
/// contig, so a mismatched count has to be a miss rather than a silent shift.
#[test]
fn a_table_of_the_wrong_length_is_refused() {
    let home = tempfile::tempdir().unwrap();
    let path = home.path().join("families.gz");
    write(&path, &annotation()).unwrap();

    assert!(read(&path, 3).is_err());
}

/// The key has to move when any of the three things the tables depend on moves, or a sweep
/// reuses one dataset's gene families for another.
#[test]
fn the_key_follows_the_contigs_and_the_database() {
    let home = std::path::Path::new("/cache");
    let names = ["one".to_string(), "two".to_string()];
    let database = std::path::Path::new("uniref100.dmnd");

    let held = path_for(home, "assembly.fa", &names, database);
    assert_eq!(held, path_for(home, "assembly.fa", &names, database));
    assert_ne!(held, path_for(home, "other.fa", &names, database));
    assert_ne!(
        held,
        path_for(home, "assembly.fa", &names[..1], database)
    );
    assert_ne!(
        held,
        path_for(home, "assembly.fa", &names, std::path::Path::new("other.dmnd"))
    );
}
