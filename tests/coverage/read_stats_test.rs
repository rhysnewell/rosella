use anyhow::Result;
use rosella::coverage::depth_file::{self, Contigs, Depth};
use rosella::coverage::read_stats::{read_counts, tie_share, tied_reads};

fn write_ties(directory: &std::path::Path, sample: &str, triples: &[[u32; 3]]) -> Result<()> {
    let bytes = triples
        .iter()
        .flatten()
        .flat_map(|word| word.to_le_bytes())
        .collect::<Vec<_>>();
    let path = directory.join(format!("{sample}.ties.zst"));
    std::fs::write(path, zstd::encode_all(bytes.as_slice(), 3)?)?;
    Ok(())
}

#[test]
fn tied_reads_sum_pairs_over_samples_for_wanted_contigs_only() -> Result<()> {
    let directory = tempfile::tempdir()?;
    write_ties(directory.path(), "a", &[[0, 1, 5], [1, 3, 2]])?;
    write_ties(directory.path(), "b", &[[0, 2, 1]])?;
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

// The largest partner, not the sum over partners, stands for a contig's tied reads, and rows
// come back in the order the names ask for, one count per sample.
#[test]
fn read_counts_take_half_the_largest_partner_off_when_untied() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let contigs = Contigs {
        names: ["a", "b", "c"].map(String::from).to_vec(),
        lengths: vec![1000; 3],
    };
    contigs.keep(directory.path())?;
    for (sample, reads) in [("s1", vec![10, 4, 7]), ("s2", vec![3, 0, 5])] {
        let depth = Depth {
            means: vec![0.0; 3],
            variances: vec![0.0; 3],
            reads,
        };
        depth_file::write(&depth_file::path(directory.path(), sample), &depth, &[])?;
    }
    write_ties(directory.path(), "s1", &[[0, 1, 6], [0, 2, 2]])?;
    let samples = ["s1", "s2"].map(String::from);
    let names = ["c", "a"].map(String::from);
    let untied = read_counts(directory.path(), &samples, &names, true)?;
    assert_eq!(untied, [6.0, 5.0, 7.0, 3.0]);
    let placed = read_counts(directory.path(), &samples, &names, false)?;
    assert_eq!(placed, [7.0, 5.0, 10.0, 3.0]);
    Ok(())
}
