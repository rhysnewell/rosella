use std::collections::{BTreeMap, HashSet};

use ndarray::array;

use rosella::embedding::knn::KnnGraph;
use rosella::refine::strays::{Tie, adrift, strays};

#[test]
fn a_tied_twin_at_another_depth_leaves_and_a_repeat_at_the_bin_depth_stays() {
    let lengths = [20_000, 20_000, 20_000, 3_000, 3_000, 3_000, 3_000];
    let depths = array![
        [10.0, 2.0],
        [11.0, 2.0],
        [10.0, 2.2],
        [10.0, 2.0],
        [40.0, 9.0],
        [10.5, 2.1],
        [9.5, 1.9],
    ];
    let bins = [vec![0, 1, 2, 3, 4, 5, 6]];
    let ties = [
        Tie {
            first: 3,
            second: 4,
            reads: 400,
        },
        Tie {
            first: 5,
            second: 6,
            reads: 400,
        },
        Tie {
            first: 0,
            second: 4,
            reads: 400,
        },
    ];
    let reads = [5000, 5000, 5000, 600, 600, 600, 600];
    let out = strays(
        bins.iter().map(Vec::as_slice),
        &lengths,
        depths.view(),
        &ties,
        &reads,
    );
    assert_eq!(out, vec![4]);
}

#[test]
fn a_short_contig_leaves_an_unfinished_bin_whose_neighbours_mostly_sit_elsewhere() {
    let lengths = [20_000, 20_000, 20_000, 20_000, 3_000, 3_000, 20_000];
    let bins = BTreeMap::from([(0, vec![0, 1, 4, 5, 6]), (1, vec![2, 3])]);
    let knn = KnnGraph {
        indices: array![[1, 4], [0, 4], [3, 4], [2, 4], [2, 3], [0, 1], [2, 3]],
        dists: array![
            [0.1, 0.1],
            [0.1, 0.1],
            [0.1, 0.1],
            [0.1, 0.1],
            [0.1, 0.1],
            [0.1, 0.1],
            [0.1, 0.1]
        ],
    };
    assert_eq!(adrift(&bins, &HashSet::new(), &lengths, &knn), vec![4]);
    assert!(adrift(&bins, &HashSet::from([0]), &lengths, &knn).is_empty());
}
