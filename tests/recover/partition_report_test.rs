use std::collections::BTreeMap;

use rosella::clustering::clusterer::Partitioning;
use rosella::clustering::graph_partition::Partition;
use rosella::recover::partition_report::{PartitionReport, keep_combine, pass_file};

fn partitioning(groups: &[&[usize]]) -> Partitioning {
    Partitioning {
        cluster_map: groups
            .iter()
            .enumerate()
            .map(|(label, members)| (label, members.to_vec()))
            .collect::<BTreeMap<_, _>>(),
        outliers: Vec::new(),
        score: None,
        arm: Partition::Leiden,
        seed: 0,
    }
}

/// Readers take the derived block as the partition the run used, so a refused pass written
/// after the kept one must not carry that name.
#[test]
fn only_the_kept_pass_reads_as_derived() {
    let home = tempfile::tempdir().unwrap();
    let mut report = PartitionReport::new(&[0, 1]);
    for (pass, groups) in [("pass1", [&[0, 1][..]]), ("pass2", [&[0][..]])] {
        report.pass(pass);
        report.add("combined", &partitioning(&groups));
    }
    report.keep("pass1");
    let written = home.path().join("partitions.tsv");
    let names = ["a".to_string(), "b".to_string()];
    report.write(&written, &names, &[10, 20]).unwrap();
    let text = std::fs::read_to_string(&written).unwrap();
    let mut lines = text.lines();
    assert_eq!(
        lines.next(),
        Some("contig\tlength\tderived_combined\tpass2_combined")
    );
    assert_eq!(lines.next(), Some("a\t10\t0\t0"));
    assert_eq!(lines.next(), Some("b\t20\t0\t-"));

    let combine = home.path().join("combine.tsv");
    let passes = ["sample".to_string(), "pass1".to_string()];
    for pass in &passes {
        std::fs::write(pass_file(&combine, pass), pass).unwrap();
    }
    keep_combine(&combine, "pass1", &passes).unwrap();
    assert_eq!(std::fs::read_to_string(&combine).unwrap(), "pass1");
    assert!(!pass_file(&combine, "sample").exists());
}
