use std::collections::HashSet;

use anyhow::Result;
use rosella::coverage::read_stats::tied_pairs;

#[test]
fn tied_pairs_join_samples_keep_wanted_contigs_and_skip_missing_files() -> Result<()> {
    let directory = tempfile::tempdir()?;
    for (sample, triples) in [("a", vec![[0u32, 1, 5], [1, 3, 2]]), ("b", vec![[0, 2, 1]])] {
        let bytes = triples
            .iter()
            .flatten()
            .flat_map(|word| word.to_le_bytes())
            .collect::<Vec<_>>();
        let path = directory.path().join(format!("{sample}.ties.zst"));
        std::fs::write(path, zstd::encode_all(bytes.as_slice(), 3)?)?;
    }
    let samples = ["a", "b", "c"].map(String::from);
    let wanted = [true, true, true, false];
    let pairs = tied_pairs(directory.path(), &samples, &wanted)?;
    assert_eq!(pairs, HashSet::from([(0, 1), (0, 2)]));
    Ok(())
}
