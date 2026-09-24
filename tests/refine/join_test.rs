use std::collections::BTreeMap;

use ndarray::Array2;
use rosella::assembly_graph::Link;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::join::{JoinSettings, join};

use crate::scorer::GenomeScorer;

#[test]
fn a_graph_lets_join_fuse_only_the_halves_a_link_joins() {
    let lengths = vec![100_000; 4];
    let coverage = Array2::zeros((lengths.len(), 2));
    let tnf = Array2::zeros((lengths.len(), 2));
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let scorer = GenomeScorer::new(vec![0, 0, 1, 1]);
    let settings = JoinSettings {
        completeness: 90.0,
        contamination: 5.0,
        max_bin_size: 15_000_000,
    };
    let link = |from, to| Link {
        from,
        to,
        branching: 1,
        walked: false,
    };
    let halves = || BTreeMap::from([(0, vec![0]), (1, vec![1]), (2, vec![2]), (3, vec![3])]);

    for (links, joined) in [
        (None, 2),
        (Some(vec![]), 0),
        (Some(vec![link(0, 1), link(1, 2)]), 1),
    ] {
        let mut bins = halves();
        let ledger = join(&features, &scorer, &mut bins, settings, links.as_deref());
        assert_eq!(ledger.joined, joined, "{links:?}");
    }
}
