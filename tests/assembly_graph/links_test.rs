//! The file is mostly sequence and the depth tables have already dropped short contigs, so the
//! reader has to skip both without the caller noticing.

use rosella::assembly_graph::{Link, read_links};

fn names(count: usize) -> Vec<String> {
    (1..=count).map(|at| format!("edge_{at}")).collect()
}

fn pairs(links: &[Link]) -> Vec<(usize, usize)> {
    links.iter().map(|link| (link.from, link.to)).collect()
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
fn a_gzipped_graph_reads_the_same_links() {
    let plain = std::fs::read("tests/data/links_branching.gfa").unwrap();
    let packed = tempfile::Builder::new()
        .suffix(".gfa.gz")
        .tempfile()
        .unwrap();
    let mut encoder = flate2::write::GzEncoder::new(packed.as_file(), flate2::Compression::fast());
    std::io::Write::write_all(&mut encoder, &plain).unwrap();
    encoder.finish().unwrap();

    assert_eq!(
        read_links(packed.path(), &names(15)).unwrap(),
        read_links("tests/data/links_branching.gfa", &names(15)).unwrap()
    );
}
