use std::path::Path;

use rosella::external::hmmer_engine::search_held;

const FLOOR: f64 = 6.66;

fn rows(table: &str) -> Vec<Vec<&str>> {
    let mut rows = table
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split_whitespace().collect::<Vec<_>>())
        .collect::<Vec<_>>();
    rows.sort();
    rows
}

// Each table was written by the hmmsearch 3.4 binary on that CPU family, so any drift in how
// proteins are held, sliced or scored shows up as a changed row. x86 flushes tiny E-values to 0.
#[test]
fn the_held_search_reports_what_hmmsearch_reports_in_one_slice_or_many() {
    rosella::pool::init(2).ok();
    let directory = tempfile::tempdir().unwrap();
    let hmm = directory.path().join("markers.hmm");
    let mut packed =
        flate2::read::GzDecoder::new(std::fs::File::open("data/gtdb_markers.hmm.gz").unwrap());
    std::io::copy(&mut packed, &mut std::fs::File::create(&hmm).unwrap()).unwrap();
    let proteins = Path::new("tests/data/hmmer/proteins.faa");
    let expected = std::fs::read_to_string(format!(
        "tests/data/hmmer/expected_{}.tbl",
        std::env::consts::ARCH
    ))
    .unwrap();
    for held in [usize::MAX, 8_000] {
        let table = search_held(&hmm, proteins, FLOOR, held).unwrap();
        assert_eq!(rows(&table), rows(&expected), "held {held} residues");
    }
}
