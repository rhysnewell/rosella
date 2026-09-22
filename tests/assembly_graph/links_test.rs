//! The file is mostly sequence and the depth tables have already dropped short contigs, so the
//! reader has to skip both without the caller noticing.

use rosella::assembly_graph::{Link, read_links};

fn names(count: usize) -> Vec<String> {
    (1..=count).map(|at| format!("edge_{at}")).collect()
}

fn pairs(links: &[Link]) -> Vec<(usize, usize)> {
    links.iter().map(|link| (link.from, link.to)).collect()
}

fn found(pair: (usize, usize), links: &[Link]) -> Link {
    *links
        .iter()
        .find(|link| (link.from, link.to) == pair)
        .expect("the pair is in the graph")
}

#[test]
fn links_are_deduped_undirected_pairs_of_surviving_contigs() {
    let links = read_links("tests/data/links.gfa", &names(3)).expect("the fixture reads");
    assert_eq!(
        pairs(&links),
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
fn branching_counts_the_busier_of_the_two_ends() {
    let links = read_links("tests/data/links_branching.gfa", &names(15)).expect("the fixture reads");

    assert_eq!(found((0, 1), &links).branching, 1, "neither end branches");
    assert_eq!(
        found((10, 12), &links).branching,
        4,
        "one end of the hub offers four continuations"
    );
}

#[test]
fn a_contig_path_marks_the_pair_it_crossed_and_no_other() {
    let links = read_links("tests/data/links_branching.gfa", &names(15)).expect("the fixture reads");

    assert!(found((10, 11), &links).walked);
    assert!(
        !found((10, 12), &links).walked,
        "a sibling off the same hub"
    );
}
