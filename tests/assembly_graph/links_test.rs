//! The file is mostly sequence and the depth tables have already dropped short contigs, so the
//! reader has to skip both without the caller noticing.

use rosella::assembly_graph::{Link, read_links};

fn names(count: usize) -> Vec<String> {
    (1..=count).map(|at| format!("edge_{at}")).collect()
}

fn pairs(links: &[Link]) -> Vec<(usize, usize)> {
    links.iter().map(|link| (link.from, link.to)).collect()
}

fn trust(links: &[Link], pair: (usize, usize)) -> f32 {
    links
        .iter()
        .find(|link| (link.from, link.to) == pair)
        .map(|link| link.trust)
        .expect("the pair is in the graph")
}

#[test]
fn links_are_deduped_undirected_pairs_of_surviving_contigs() {
    let found = read_links("tests/data/links.gfa", &names(3)).expect("the fixture reads");
    assert_eq!(
        pairs(&found),
        vec![(0, 1), (0, 2)],
        "the reciprocal pair collapses to one, the self link and the link to a filtered contig go"
    );
}

#[test]
fn a_contig_the_filter_dropped_takes_its_links_with_it() {
    let held = vec!["edge_2".to_string(), "edge_3".to_string()];
    assert!(
        read_links("tests/data/links.gfa", &held)
            .expect("the fixture reads")
            .is_empty()
    );
}

#[test]
fn a_branching_end_is_trusted_below_an_unbranched_one() {
    let found = read_links("tests/data/links_trust.gfa", &names(15)).expect("the fixture reads");

    assert_eq!(
        trust(&found, (0, 1)),
        1.0,
        "both ends offer one continuation"
    );
    assert_eq!(
        trust(&found, (10, 12)),
        0.25,
        "one end of the hub offers four, and the median link in this graph offers one"
    );
}

#[test]
fn a_contig_path_floors_a_branching_link_at_the_median() {
    let found = read_links("tests/data/links_trust.gfa", &names(15)).expect("the fixture reads");

    assert_eq!(
        trust(&found, (10, 11)),
        1.0,
        "the assembler walked a contig through this pair, so the branching does not discount it"
    );
}
