use ndarray::array;
use rosella::coverage::coverage_table::CoverageTable;
use rosella::recover::abundance;

#[test]
fn depth_is_weighted_by_length_and_shares_split_the_binned_depth() {
    let coverage = CoverageTable {
        table: array![[2.0, 0.1], [4.0, 0.1], [1.0, 0.1]],
        average_depths: vec![2.0, 4.0, 1.0],
        contig_names: vec!["a".into(), "b".into(), "c".into()],
        contig_lengths: vec![100, 300, 600],
        sample_names: vec!["reads".into()],
    };
    let bins = vec![
        ("one".to_string(), vec![0, 1]),
        ("two".to_string(), vec![2]),
    ];
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join(abundance::ABUNDANCE_FILE);

    abundance::write(&bins, &coverage, &path).unwrap();

    let written = std::fs::read_to_string(path).unwrap();
    let rows = written.lines().collect::<Vec<_>>();
    assert_eq!(rows[0], "bin\tbp\treads_depth\treads_share");
    assert_eq!(rows[1], "one\t400\t3.5000\t70.0000");
    assert_eq!(rows[2], "two\t600\t1.0000\t30.0000");
}
