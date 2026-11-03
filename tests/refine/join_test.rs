use std::collections::BTreeMap;

use ndarray::{Array2, s};
use rosella::embedding::features::ContigFeatures;
use rosella::refine::join::{Depths, JoinSettings, join};

use crate::scorer::MarkerScorer;

const SETTINGS: JoinSettings = JoinSettings {
    completeness: 90.0,
    contamination: 5.0,
    max_bin_size: 10_000_000,
};

// Both pieces complete the host's markers, and markers alone take the first, off-depth one.
#[test]
fn a_piece_joins_only_the_host_whose_depth_it_shares() {
    let jitter = [1.0, 1.1, 0.9, 1.05, 0.95];
    let profiles = [[10.0, 20.0, 5.0], [10.0, 2.0, 40.0], [10.0, 20.0, 5.0]];
    let mut depths = Array2::zeros((15, 3));
    for (group, profile) in profiles.iter().enumerate() {
        for (at, scale) in jitter.iter().enumerate() {
            for (sample, depth) in profile.iter().enumerate() {
                depths[[group * 5 + at, sample]] = depth * scale;
            }
        }
    }
    let piece = [vec![5], vec![6], vec![7], vec![8], vec![]];
    let markers = (0..5)
        .map(|marker| vec![marker])
        .chain(piece.iter().cloned())
        .chain(piece.iter().cloned())
        .collect::<Vec<_>>();
    let quality = MarkerScorer::new(markers, 9);
    let lengths = vec![10_000; 15];
    let (coverage, tnf) = (Array2::zeros((15, 2)), Array2::zeros((15, 2)));
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let start = (0..3)
        .map(|group| (group, (group * 5..group * 5 + 5).collect::<Vec<_>>()))
        .collect::<BTreeMap<_, _>>();

    let joined = |table| {
        let mut bins = start.clone();
        join(
            &features,
            &quality,
            Depths {
                table,
                lengths: &lengths,
            },
            &mut bins,
            SETTINGS,
        );
        bins
    };
    let all = joined(depths.view());
    assert_eq!(all.len(), 2);
    assert_eq!(all[&0], (0..5).chain(10..15).collect::<Vec<_>>());
    assert_eq!(
        joined(depths.slice(s![.., ..1])),
        start,
        "one sample cannot tell the pieces apart"
    );
}
