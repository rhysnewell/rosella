use ndarray::{Array2, s};
use rosella::recover::coherence::{Coherence, tie_spread};

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
            .holds_one(
                members,
                depths,
                &lengths,
                &[150.0; 3][..depths.ncols()],
                |_, _| 0.0,
            )
    };
    assert!(holds(depths.view(), &one));
    assert!(!holds(depths.view(), &mixed));
    assert!(
        holds(depths.slice(s![.., ..1]), &mixed),
        "both genomes sit at one depth in the first sample"
    );
}

#[test]
fn a_low_depth_genome_holds_once_its_reads_are_counted() {
    let jitter = [1.0, 1.1, 0.9, 1.05, 0.95];
    let mut depths = Array2::zeros((21, 2));
    for (group, profile) in [[10.0, 20.0], [30.0, 2.0], [4.0, 9.0]].iter().enumerate() {
        for (at, scale) in jitter.iter().enumerate() {
            for (sample, depth) in profile.iter().enumerate() {
                depths[[group * 5 + at, sample]] = depth * scale;
            }
        }
    }
    for (at, sparse) in [0.0, 0.5, 0.0, 0.5, 0.5, 0.0].iter().enumerate() {
        depths[[15 + at, 0]] = 1.0;
        depths[[15 + at, 1]] = *sparse;
    }
    let hosts = (0..3)
        .map(|group| (group * 5..group * 5 + 5).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let hosts = hosts.iter().map(Vec::as_slice).collect::<Vec<_>>();
    let lengths = vec![5_000; 21];
    let coherence = Coherence::of(&hosts, depths.view(), &lengths).expect("three hosts");
    let sparse = (15..21).collect::<Vec<_>>();
    assert!(coherence.holds_one(
        &sparse,
        depths.view(),
        &lengths,
        &[4_000.0, 4_000.0],
        |_, _| 0.0
    ));
    assert!(
        !coherence.holds_one(&sparse, depths.view(), &lengths, &[150.0, 150.0], |_, _| {
            0.0
        })
    );
}

#[test]
fn a_contig_the_tie_split_starved_holds_once_its_tie_spread_counts() {
    let jitter = [1.0, 1.1, 0.9, 1.05, 0.95];
    let mut depths = Array2::zeros((19, 3));
    let profiles = [[10.0, 20.0, 5.0], [30.0, 2.0, 40.0], [4.0, 9.0, 60.0]];
    for (group, profile) in profiles.iter().enumerate() {
        for (at, scale) in jitter.iter().enumerate() {
            for (sample, depth) in profile.iter().enumerate() {
                depths[[group * 5 + at, sample]] = depth * scale;
            }
        }
    }
    for (at, scale) in [1.0, 1.05, 0.95, 1.0].iter().enumerate() {
        for (sample, depth) in [10.0, 80.0, 1.0].iter().enumerate() {
            depths[[15 + at, sample]] = depth * scale;
        }
    }
    depths[[18, 1]] = 20.0;
    let hosts = (0..3)
        .map(|group| (group * 5..group * 5 + 5).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let hosts = hosts.iter().map(Vec::as_slice).collect::<Vec<_>>();
    let lengths = vec![10_000; 19];
    let coherence = Coherence::of(&hosts, depths.view(), &lengths).expect("three hosts");
    let bin = [15, 16, 17, 18];
    let spread = tie_spread(4_000, 5_333, 150.0, 10_000);
    let tied = |contig, sample| match (contig, sample) {
        (18, 1) => spread,
        _ => 0.0,
    };
    assert!(!coherence.holds_one(&bin, depths.view(), &lengths, &[150.0; 3], |_, _| 0.0));
    assert!(coherence.holds_one(&bin, depths.view(), &lengths, &[150.0; 3], tied));
}
