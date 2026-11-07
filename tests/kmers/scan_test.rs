use std::fs;
use std::io::Write;
use std::path::Path;

use rosella::kmers::scan::{Floors, scan};

fn content(path: &Path) -> u64 {
    let floors = Floors {
        halves: Some(usize::MAX),
        ..Floors::default()
    };
    scan(path.to_str().unwrap(), 4, floors)
        .unwrap()
        .content
        .unwrap()
}

#[test]
fn a_reordered_wrapped_or_packed_copy_reads_as_the_same_contigs() {
    let home = tempfile::tempdir().unwrap();
    let write = |file: &str, text: &str| {
        let path = home.path().join(file);
        fs::write(&path, text).unwrap();
        path
    };
    let original = write("original.fasta", ">a one\nACGTACGTAA\n>b\nGGGCCC\n");
    let reordered = write(
        "reordered.fasta",
        ">b two\nGGG\r\nCCC\r\n>a\nACGTA\nCGTAA\n",
    );
    let packed = home.path().join("packed.fasta.gz");
    let mut encoder = flate2::write::GzEncoder::new(
        fs::File::create(&packed).unwrap(),
        flate2::Compression::fast(),
    );
    encoder.write_all(&fs::read(&original).unwrap()).unwrap();
    encoder.finish().unwrap();

    let wanted = content(&original);
    assert_eq!(
        wanted, 16520860808620239887,
        "a changed digest renames every cache entry"
    );
    assert_eq!(content(&reordered), wanted);
    assert_eq!(content(&packed), wanted);
    for other in [
        ">a\nACGTACGTAT\n>b\nGGGCCC\n",
        ">a\nACGTACGTAA\n>c\nGGGCCC\n",
        ">a\nACGTACGTAA\n",
    ] {
        assert_ne!(content(&write("other.fasta", other)), wanted, "{other}");
    }
}
