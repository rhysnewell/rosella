//! The sweep bound was fixed at 11 for every assembly. Raising it turned out to change no
//! result, because DBCV discards the larger values, so what is left to pin is the cost: the
//! bound is reachable and moving it must not multiply the number of HDBSCAN fits.

use rosella::clustering::clusterer::{DEFAULT_LARGEST_CLUSTER, SWEEP_WIDTH, cluster_sizes};

#[test]
fn the_default_bound_is_the_run_of_integers_flight_swept() {
    assert_eq!(
        cluster_sizes(DEFAULT_LARGEST_CLUSTER),
        (2..=DEFAULT_LARGEST_CLUSTER).collect::<Vec<_>>()
    );
}

#[test]
fn a_raised_bound_still_costs_ten_points() {
    for largest in [50, 200, 2000] {
        let sizes = cluster_sizes(largest);
        assert!(
            sizes.len() <= SWEEP_WIDTH,
            "{largest} gave {} points",
            sizes.len()
        );
        assert_eq!(sizes.first().copied(), Some(2));
        assert_eq!(sizes.last().copied(), Some(largest));
        assert!(sizes.windows(2).all(|pair| pair[0] < pair[1]), "{sizes:?}");
    }
}

/// A bound below the smallest cluster HDBSCAN can form would leave nothing to sweep.
#[test]
fn a_bound_under_two_still_sweeps_something() {
    assert_eq!(cluster_sizes(0), vec![2]);
}
