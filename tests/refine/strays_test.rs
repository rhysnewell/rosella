use std::collections::BTreeMap;

use rosella::refine::strays::gather;

use crate::scorer::MarkerScorer;

const LENGTHS: [usize; 5] = [5_000, 4_000, 1_000, 3_000, 800];

#[test]
fn strays_on_a_bridge_chain_go_to_the_bin_holding_most_of_it() {
    let cases = [
        (
            "a stray and an unbinned contig move, the stray's bin keeps the rest",
            vec![vec![0], vec![1], vec![2], vec![3], vec![4]],
            vec![(1, 2, 1), (1, 4, 3)],
            BTreeMap::from([(0, vec![0, 1, 2, 4]), (1, vec![3])]),
            vec![],
        ),
        (
            "a stray that repeats a home marker stays put",
            vec![vec![0], vec![1], vec![1], vec![3], vec![4]],
            vec![(1, 2, 1)],
            BTreeMap::from([(0, vec![0, 1]), (1, vec![2, 3])]),
            vec![4],
        ),
        (
            "a bin left empty is dropped",
            vec![vec![0], vec![1], vec![2], vec![3], vec![4]],
            vec![(1, 2, 1), (1, 3, 1)],
            BTreeMap::from([(0, vec![0, 1, 2, 3])]),
            vec![4],
        ),
    ];
    for (case, markers, bridges, bins, unbinned) in cases {
        let scorer = MarkerScorer::new(markers, 5);
        let mut held = BTreeMap::from([(0, vec![0, 1]), (1, vec![2, 3])]);
        let mut loose = vec![4];
        gather(&mut held, &mut loose, &bridges, &LENGTHS, &scorer);
        assert_eq!((held, loose), (bins, unbinned), "{case}");
    }
}
