use std::io::Write;

use rosella::quality::orfs::call_over;

fn assembly(contigs: &[(&str, usize)]) -> tempfile::NamedTempFile {
    let mut file = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
    let unit = b"ATGAAACGTTTAGGCTCAGAAGTTTGCGACGCTATCCGTGAATTCTAA";
    for (name, repeats) in contigs {
        writeln!(file, ">{name}").unwrap();
        for _ in 0..*repeats {
            file.write_all(unit).unwrap();
        }
        writeln!(file).unwrap();
    }
    file.flush().unwrap();
    file
}

#[test]
fn orfs_arrive_indexed_by_the_contig_they_came_from() {
    rosella::pool::init(2).ok();
    let file = assembly(&[("one", 80), ("two", 120), ("three", 60)]);
    let mut seen: Vec<usize> = Vec::new();
    let walked = call_over(file.path().to_str().unwrap(), 100..usize::MAX, |batch| {
        seen.extend(batch.iter().map(|orf| orf.contig));
        Ok(())
    })
    .unwrap();
    assert_eq!(walked.names, ["one", "two", "three"]);
    assert!(seen.iter().all(|contig| *contig < walked.names.len()));
    let mut ordered = seen.clone();
    ordered.sort_unstable();
    assert_eq!(seen, ordered, "contigs are reported in assembly order");
}

#[test]
fn a_contig_past_the_band_is_named_but_never_called() {
    rosella::pool::init(2).ok();
    let file = assembly(&[("one", 80), ("two", 120), ("three", 60)]);
    let mut seen: Vec<usize> = Vec::new();
    let walked = call_over(file.path().to_str().unwrap(), 100..4000, |batch| {
        seen.extend(batch.iter().map(|orf| orf.contig));
        Ok(())
    })
    .unwrap();
    assert_eq!(walked.names, ["one", "two", "three"]);
    assert_eq!(walked.lengths, [3840, 5760, 2880]);
    seen.dedup();
    assert_eq!(
        seen,
        [0, 2],
        "the long contig keeps its slot and is skipped"
    );
}

#[test]
fn an_error_from_the_callback_is_what_comes_back() {
    rosella::pool::init(2).ok();
    let file = assembly(&[("one", 80)]);
    let failed = call_over(file.path().to_str().unwrap(), 100..usize::MAX, |_| {
        Err(anyhow::anyhow!("the sink said no"))
    });
    assert_eq!(failed.unwrap_err().to_string(), "the sink said no");
}

#[test]
fn a_missing_assembly_surfaces_the_readers_error() {
    rosella::pool::init(2).ok();
    let failed = call_over("/nowhere/at/all.fa", 100..usize::MAX, |_| Ok(()));
    assert!(failed.is_err());
}
