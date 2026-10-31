use std::collections::{BTreeMap, HashSet};

use ndarray::array;

use rosella::embedding::knn::KnnGraph;
use rosella::refine::adrift::adrift;

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
