use std::collections::BTreeMap;

use rosella::clustering::clusterer::Partitioning;
use rosella::clustering::graph_partition::Partition;
use rosella::recover::partition_report::{PartitionReport, keep_pass, pass_file, write_pull};

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
    keep_pass(&combine, "pass1", &passes).unwrap();
    assert_eq!(std::fs::read_to_string(&combine).unwrap(), "pass1");
    assert!(!pass_file(&combine, "sample").exists());
}

/// Rows are graph positions and the names are contigs, so a mixup between the two misnames every
/// row, and a neighbour outside every community must still be counted.
#[test]
fn pull_sums_each_community_with_and_without_lengths() {
    let home = tempfile::tempdir().unwrap();
    let graph = sprs::CsMatI::<f32, u32, usize>::new(
        (3, 3),
        vec![0, 2, 3, 4],
        vec![1, 2, 0, 0],
        vec![0.5, 0.25, 0.5, 0.25],
    );
    let mut held = partitioning(&[&[0, 1]]);
    held.outliers = vec![2];
    let names = ["x", "a", "b", "c"].map(str::to_string);
    let written = home.path().join("pull.tsv");
    write_pull(
        &written,
        &graph,
        &held,
        &[1, 2, 3],
        &names,
        &[0, 100, 400, 900],
    )
    .unwrap();
    let text = std::fs::read_to_string(&written).unwrap();
    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        [
            "contig\tcommunity\ttoward\tweight\tbp_weight",
            "a\t0\t-\t0.250000\t75.0",
            "a\t0\t0\t0.500000\t100.0",
            "b\t0\t0\t0.500000\t100.0",
            "c\t-\t0\t0.250000\t75.0",
        ]
    );
}
