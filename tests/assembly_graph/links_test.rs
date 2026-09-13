//! The file is mostly sequence and the depth tables have already dropped short contigs, so the
//! reader has to skip both without the caller noticing.

use rosella::assembly_graph::read_links;

fn names() -> Vec<String> {
    ["edge_1", "edge_2", "edge_3"]
        .iter()
        .map(|name| name.to_string())
        .collect()
}

#[test]
fn links_are_deduped_undirected_pairs_of_surviving_contigs() {
    let found = read_links("tests/data/links.gfa", &names()).expect("the fixture reads");
    assert_eq!(
        found,
        vec![(0, 1), (0, 2)],
        "the reciprocal pair collapses to one, the self link and the link to a filtered contig go"
    );
}

#[test]
fn a_contig_the_filter_dropped_takes_its_links_with_it() {
    let held = vec!["edge_2".to_string(), "edge_3".to_string()];
    assert!(read_links("tests/data/links.gfa", &held).expect("the fixture reads").is_empty());
}
