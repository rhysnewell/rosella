use std::io::{BufRead, BufReader};
use std::path::Path;

use rosella::external::hmmer_engine::{align, search_held};

const FLOOR: f64 = 6.66;

fn markers() -> impl Iterator<Item = String> {
    let packed = std::fs::File::open("data/gtdb_markers.hmm.gz").unwrap();
    BufReader::new(flate2::read::GzDecoder::new(packed))
        .lines()
        .map(Result::unwrap)
}

fn profile(name: &str) -> String {
    let mut block = String::new();
    for line in markers() {
        block += &line;
        block.push('\n');
        if line.starts_with("//") {
            if block
                .lines()
                .any(|held| held.split_whitespace().eq(["NAME", name]))
            {
                return block;
            }
            block.clear();
        }
    }
    panic!("no {name} model");
}

fn rows(table: &str) -> Vec<Vec<&str>> {
    let mut rows = table
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split_whitespace().collect::<Vec<_>>())
        .collect::<Vec<_>>();
    rows.sort();
    rows
}

// The table was written by the hmmsearch 3.4 binary on Linux x86_64, which flushes tiny E-values
// to 0. The search computes the same bits on every CPU, so one table serves both.
#[test]
fn the_held_search_reports_what_hmmsearch_reports_in_one_slice_or_many() {
    rosella::pool::init(2).ok();
    let directory = tempfile::tempdir().unwrap();
    let hmm = directory.path().join("markers.hmm");
    std::fs::write(&hmm, markers().map(|line| line + "\n").collect::<String>()).unwrap();
    let proteins = Path::new("tests/data/hmmer/proteins.faa");
    let expected = std::fs::read_to_string("tests/data/hmmer/expected.tbl").unwrap();
    for held in [usize::MAX, 8_000] {
        let table = search_held(&hmm, proteins, FLOOR, held).unwrap();
        assert_eq!(rows(&table), rows(&expected), "held {held} residues");
    }
}

// The expected rows were written by the hmmalign 3.4 binary, alike on ARM and x86.
#[test]
fn the_held_alignment_writes_what_hmmalign_writes() {
    let fasta = std::fs::read_to_string("tests/data/hmmer/era.faa").unwrap();
    let proteins = fasta.lines().skip(1).step_by(2).collect::<Vec<_>>();
    let expected = std::fs::read_to_string("tests/data/hmmer/era.a2m").unwrap();
    assert_eq!(align(&profile("era"), &proteins).unwrap(), expected);
}
