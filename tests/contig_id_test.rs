use rosella::contig_id;

#[test]
fn an_annotated_header_matches_the_bare_depth_table_name() {
    assert_eq!(contig_id(b"edge_1 LN:i:29210 RC:i:1838").unwrap(), "edge_1");
    assert_eq!(contig_id(b"edge_1").unwrap(), "edge_1");
    assert_eq!(contig_id(b"edge_1\tlength=900").unwrap(), "edge_1");
    assert_eq!(contig_id(b"").unwrap(), "");
}
