use anyhow::Result;
use rosella::coverage::read_stats::{tie_share, tied_reads};

#[test]
fn tied_reads_sum_pairs_over_samples_for_wanted_contigs_only() -> Result<()> {
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
    let tied = tied_reads(directory.path(), &samples, &[true, true, false, false])?;
    assert_eq!(tied, [6, 7, 0, 0]);
    Ok(())
}

// Two contigs that share every read, primaries split evenly, and a contig with half its reads
// tied to one other contig.
#[test]
fn tie_share_is_exact_for_reads_tied_between_two_contigs() {
    assert_eq!(tie_share(100, 50), 1.0);
    assert_eq!(tie_share(50, 75), 0.5);
    assert_eq!(tie_share(0, 0), 0.0);
}
