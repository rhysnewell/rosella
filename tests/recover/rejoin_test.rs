use ndarray::Array2;
use rosella::embedding::knn::KnnGraph;
use rosella::recover::fold::{DepthFit, Role};
use rosella::recover::rejoin::rejoins;

#[test]
fn pieces_rejoin_their_largest_sibling_at_their_depth_only() {
    let ranges = [0..5, 5..10, 10..16, 16..21, 21..24];
    let profiles = [
        [10.0, 20.0],
        [5.0, 8.0],
        [5.0, 8.0],
        [30.0, 2.0],
        [5.0, 8.0],
    ];
    let jitter = [1.0, 1.1, 0.9, 1.05, 0.95, 1.02];
    let mut depths = Array2::zeros((24, 2));
    for (range, profile) in ranges.iter().zip(profiles) {
        for (at, contig) in range.clone().enumerate() {
            for (sample, depth) in profile.iter().enumerate() {
                depths[[contig, sample]] = depth * jitter[at];
            }
        }
    }
    let near = |own: &std::ops::Range<usize>, other: &std::ops::Range<usize>, keep: usize| {
        own.clone()
            .take(keep)
            .chain(other.clone().cycle())
            .take(5)
            .map(|contig| contig as u32)
            .collect::<Vec<_>>()
    };
    let mut rows = Vec::new();
    for (group, range) in ranges.iter().enumerate() {
        let (other, keep) = match group {
            0 => (&ranges[0], 5),
            1 => (&ranges[2], 2),
            2 | 3 => (&ranges[1], 2),
            _ => (&ranges[1], 1),
        };
        for _ in range.clone() {
            rows.extend(near(range, other, keep));
        }
    }
    let knn = KnnGraph {
        indices: Array2::from_shape_vec((24, 5), rows).unwrap(),
        dists: Array2::from_elem((24, 5), 0.1),
    };
    let members = ranges
        .iter()
        .map(|range| range.clone().collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let groups = members.iter().map(Vec::as_slice).collect::<Vec<_>>();
    let roles = [
        Role::Host,
        Role::Other,
        Role::Other,
        Role::Other,
        Role::Other,
    ];
    let lengths = vec![10_000; 24];
    let depth = DepthFit::of(&groups, &roles, depths.view(), &lengths).unwrap();
    let eukaryotic = [false, true, true, true, true];
    assert_eq!(
        rejoins(&groups, &eukaryotic, &depth, &knn, &lengths, |_, _| None),
        vec![(1, 2), (4, 2)],
        "the piece at another depth stays apart, and no piece moves into one that moves"
    );
    let shared = |one: &[usize], other: &[usize]| {
        (one[0].min(other[0]) == 5 && one[0].max(other[0]) == 10).then_some(0.9)
    };
    assert_eq!(
        rejoins(&groups, &eukaryotic, &depth, &knn, &lengths, shared),
        vec![(4, 1)],
        "two pieces holding most of each other's marker families are two genomes"
    );
}
