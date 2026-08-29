//! The curve stopped being derived once the derivation was shown to be two clamps, so what
//! is left to check is which of the two ways of setting it a run picks.

use rosella::embedding::umap::{Curve, EmbedOverrides, length_weights, n_components};

#[test]
fn the_curve_is_pinned_until_min_dist_or_spread_asks_for_a_fit() {
    let pinned = |overrides: &EmbedOverrides| {
        matches!(Curve::from_overrides(overrides), Curve::Pinned(_))
    };

    assert!(pinned(&EmbedOverrides::default()));
    assert!(pinned(&EmbedOverrides {
        b: Some(0.4),
        ..Default::default()
    }));
    assert!(!pinned(&EmbedOverrides {
        spread: Some(2.0),
        ..Default::default()
    }));
    assert!(!pinned(&EmbedOverrides {
        min_dist: Some(0.1),
        ..Default::default()
    }));
}

/// A fit takes umap's own defaults for whichever of the pair was left out, rather than
/// carrying the pinned curve's values into a parameterisation they do not belong to.
#[test]
fn a_half_given_fit_completes_itself() {
    let Curve::Fit { min_dist, spread } = Curve::from_overrides(&EmbedOverrides {
        spread: Some(2.0),
        ..Default::default()
    }) else {
        panic!("spread did not ask for a fit");
    };
    assert_eq!((min_dist, spread), (0.0, 2.0));
}

#[test]
fn components_follow_the_estimate_and_fall_back_to_the_sample_count() {
    assert_eq!(n_components(Some(7.20), 1), 7);
    assert_eq!(n_components(Some(12.05), 5), 10);
    assert_eq!(n_components(Some(0.4), 4), 4);
    assert_eq!(n_components(None, 1), 2);
    assert_eq!(n_components(None, 50), 10);
}

#[test]
fn length_weights_are_empty_at_power_zero() {
    assert!(length_weights(&[1000, 2000, 4000], 0.0).is_empty());
    assert!(length_weights(&[], 1.0).is_empty());
}

/// The normalisation is what keeps the graph's maximum weight, and with it the threshold
/// that drops the weakest edges, where it was.
#[test]
fn length_weights_hold_a_geometric_mean_of_one() {
    let weights = length_weights(&[2000, 3000, 8000, 150_000], 1.0);
    let log_mean = weights.iter().map(|w| (*w as f64).ln()).sum::<f64>() / weights.len() as f64;
    assert!(
        log_mean.exp() - 1.0 < 1e-5,
        "geometric mean was {}",
        log_mean.exp()
    );
}

/// Unbounded, a megabase contig is drawn thousands of times more often than a 1.5 kb one
/// and the short contigs never move.
#[test]
fn length_weights_are_bounded_either_side() {
    let weights = length_weights(&[500, 3000, 20_000_000], 1.0);
    assert!(
        weights.iter().all(|w| (0.2..=5.0).contains(w)),
        "{weights:?}"
    );
    assert!(weights[0] < weights[1] && weights[1] < weights[2]);
}
