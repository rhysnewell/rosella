use ndarray::{Array2, s};
use rosella::recover::fold::{Role, chi_square_above, folds};

#[test]
fn chi_square_tail_reads_the_tabled_quantiles() {
    for (statistic, degrees, above) in [
        (3.841459, 1, 0.05),
        (5.991465, 2, 0.05),
        (18.307038, 10, 0.05),
        (4.865182, 10, 0.9),
    ] {
        let read = chi_square_above(statistic, degrees);
        assert!(
            (read - above).abs() < 1e-6,
            "{statistic} on {degrees}: {read}"
        );
    }
}

#[test]
fn a_donor_folds_only_where_the_samples_tell_its_host_apart() {
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
    let members = (0..3)
        .map(|group| (group * 5..group * 5 + 5).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let groups = members.iter().map(Vec::as_slice).collect::<Vec<_>>();
    let roles = [Role::Host, Role::Host, Role::Donor];
    let lengths = vec![10_000; 15];
    assert_eq!(
        folds(&groups, &roles, depths.view(), &lengths),
        vec![(2, 0)]
    );
    assert!(
        folds(&groups, &roles, depths.slice(s![.., ..1]), &lengths).is_empty(),
        "both hosts sit at the donor's depth in the first sample"
    );
}
