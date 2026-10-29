use std::path::Path;

use rand::{SeedableRng, rngs::StdRng};
use rosella::coverage::depth_file::ties_path;
use rosella::learned::partners::Partners;

fn ties(directory: &Path, sample: &str, rows: &[[u32; 3]]) {
    let bytes = rows
        .iter()
        .flatten()
        .flat_map(|word| word.to_le_bytes())
        .collect::<Vec<_>>();
    let packed = zstd::encode_all(bytes.as_slice(), 3).unwrap();
    std::fs::write(ties_path(directory, sample), packed).unwrap();
}

#[test]
fn tied_training_contigs_sit_side_by_side_and_each_once() {
    let directory = tempfile::tempdir().unwrap();
    let rows = ["a", "b", "c", "d", "e"]
        .iter()
        .map(|name| format!("{name}\t2000\n"))
        .collect::<String>();
    std::fs::write(
        directory.path().join("contigs.tsv"),
        format!("contig\tlength\n{rows}"),
    )
    .unwrap();
    ties(directory.path(), "s0", &[[0, 1, 3], [2, 3, 1]]);
    ties(directory.path(), "s1", &[[0, 1, 2], [2, 4, 9]]);
    let names = ["b", "a", "c", "e"].map(String::from);
    let samples = ["s0", "s1"].map(String::from);
    let partners = Partners::read(directory.path(), &samples, &names).unwrap();
    let mut rng = StdRng::seed_from_u64(7);
    assert_eq!(partners.paired(&[1, 2, 0, 3], &mut rng), [1, 0, 2, 3]);
    assert_eq!(partners.paired(&[3, 0, 1], &mut rng), [3, 0, 1]);
}
