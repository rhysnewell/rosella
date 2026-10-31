use std::collections::{BTreeMap, HashSet};

use ndarray::{Array2, array};

use rosella::embedding::knn::KnnGraph;
use rosella::refine::adrift::{Gate, adrift};

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
    let run = |finished| adrift(&bins, &finished, &lengths, &knn, Gate::Bars, 5.0);
    assert_eq!(run(HashSet::new()), vec![4]);
    assert!(run(HashSet::from([0])).is_empty());
}

#[test]
fn an_open_bin_sheds_only_more_adrift_contigs_than_the_finished_bins_rate_allows() {
    let drifting = [0, 1, 21]
        .into_iter()
        .chain(32..40)
        .collect::<HashSet<usize>>();
    let mut lengths = vec![3_000; 44];
    for anchor in [20, 31, 42, 43] {
        lengths[anchor] = 50_000;
    }
    let bins = BTreeMap::from([
        (0, (0..20).collect::<Vec<usize>>()),
        (1, (20..31).collect()),
        (2, (31..42).collect()),
        (3, vec![42, 43]),
    ]);
    let neighbours = (0..44u32)
        .flat_map(|contig| {
            let members = bins
                .values()
                .find(|members| members.contains(&(contig as usize)))
                .unwrap();
            match drifting.contains(&(contig as usize)) {
                true => [42, 43],
                false => {
                    let mut others = members
                        .iter()
                        .filter(|member| **member != contig as usize && !drifting.contains(member))
                        .map(|member| *member as u32);
                    let first = others.next().unwrap();
                    [first, others.next().unwrap_or(first)]
                }
            }
        })
        .collect::<Vec<_>>();
    let knn = KnnGraph {
        indices: Array2::from_shape_vec((44, 2), neighbours).unwrap(),
        dists: Array2::from_elem((44, 2), 0.1),
    };
    let finished = HashSet::from([0]);
    let shed = (32..40).collect::<Vec<usize>>();
    assert_eq!(
        adrift(&bins, &finished, &lengths, &knn, Gate::Excess, 5.0),
        shed
    );
    assert_eq!(
        adrift(&bins, &finished, &lengths, &knn, Gate::Sized, 5.0),
        shed
    );
    lengths[31] = 1_000_000;
    assert!(adrift(&bins, &finished, &lengths, &knn, Gate::Sized, 5.0).is_empty());
}
