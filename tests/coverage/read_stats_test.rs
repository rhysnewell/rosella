use anyhow::Result;
use rosella::coverage::read_stats::tie_share;

const HEADER: &str = "sample\treads\tmapped\ttied\tmean_read_length\tmapped_bases\n";

#[test]
fn tie_share_pools_samples_and_is_none_without_summaries() -> Result<()> {
    let directory = tempfile::tempdir()?;
    for (sample, mapped, tied) in [("a", 900, 90), ("b", 100, 60)] {
        let row = format!("{HEADER}{sample}\t1000\t{mapped}\t{tied}\t150.0\t1\n");
        std::fs::write(directory.path().join(format!("{sample}.summary.tsv")), row)?;
    }
    let samples = |names: &[&str]| {
        names
            .iter()
            .map(|name| name.to_string())
            .collect::<Vec<_>>()
    };
    assert_eq!(
        tie_share(directory.path(), &samples(&["a", "b"]))?,
        Some(0.15)
    );
    assert_eq!(tie_share(directory.path(), &samples(&["a", "c"]))?, None);
    Ok(())
}
