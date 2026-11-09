use ndarray::{Array2, s};
use rosella::recover::coherence::Coherence;

#[test]
fn a_bin_fails_only_where_the_samples_tell_its_genomes_apart() {
    let jitter = [1.0, 1.1, 0.9, 1.05, 0.95];
    let profiles = [
        [10.0, 20.0, 5.0],
        [30.0, 2.0, 40.0],
        [4.0, 9.0, 60.0],
        [10.0, 80.0, 1.0],
    ];
    let mut depths = Array2::zeros((20, 3));
    for (group, profile) in profiles.iter().enumerate() {
        for (at, scale) in jitter.iter().enumerate() {
            for (sample, depth) in profile.iter().enumerate() {
                depths[[group * 5 + at, sample]] = depth * scale;
            }
        }
    }
    let hosts = (0..3)
        .map(|group| (group * 5..group * 5 + 5).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let hosts = hosts.iter().map(Vec::as_slice).collect::<Vec<_>>();
    let lengths = vec![10_000; 20];
    let one = [15, 16, 17, 18];
    let mixed = [0, 1, 15, 16];
    let holds = |depths, members: &[usize]| {
        Coherence::of(&hosts, depths, &lengths)
            .expect("three hosts measure the spread")
            .holds_one(members, depths, &lengths)
    };
    assert!(holds(depths.view(), &one));
    assert!(!holds(depths.view(), &mixed));
    assert!(
        holds(depths.slice(s![.., ..1]), &mixed),
        "both genomes sit at one depth in the first sample"
    );
}
