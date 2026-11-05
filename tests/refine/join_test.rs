use std::collections::BTreeMap;

use ndarray::Array2;

use rosella::embedding::features::ContigFeatures;
use rosella::refine::join::{JoinSettings, join};

use crate::scorer::MarkerScorer;

fn bins() -> BTreeMap<usize, Vec<usize>> {
    BTreeMap::from([(0, vec![0, 1]), (1, vec![2]), (2, vec![3])])
}

#[test]
fn bridges_name_the_only_pairs_the_join_may_fuse() {
    let scorer = MarkerScorer::new(vec![vec![0], vec![1], vec![2], vec![3]], 4);
    let coverage = Array2::zeros((4, 2));
    let tnf = Array2::zeros((4, 2));
    let lengths = vec![1_000; 4];
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let settings = JoinSettings {
        completeness: 90.0,
        contamination: 10.0,
        max_bin_size: usize::MAX,
    };

    let mut bridged = bins();
    join(
        &features,
        &scorer,
        &mut bridged,
        settings,
        Some(&[(1, 3, 2)]),
    );
    assert_eq!(bridged, BTreeMap::from([(0, vec![0, 1, 3]), (1, vec![2])]));

    let mut unbridged = bins();
    join(&features, &scorer, &mut unbridged, settings, Some(&[]));
    assert_eq!(
        unbridged,
        bins(),
        "with no long reads between bins nothing joins"
    );
}
