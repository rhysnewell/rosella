use rosella::markers::genome_tree::References;

#[test]
fn nearest_reads_identity_over_shared_columns_and_breaks_ties_on_overlap() {
    rosella::pool::init(2).ok();
    let home = tempfile::tempdir().unwrap();
    let path = home.path().join("references.fasta");
    std::fs::write(&path, ">short\nAAAA--\n>long\nAAAAAA\n>far\nCCCCCC\n").unwrap();
    let references = References::read(&path).unwrap();
    assert_eq!(references.nearest(b"AAAAAC"), Some("short"));
    assert_eq!(references.nearest(b"AAAAA-"), Some("long"));
    assert_eq!(references.nearest(b"------"), None);
}
