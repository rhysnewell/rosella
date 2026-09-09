//! Half the curve stopped being derived once `b`'s derivation was shown to be a clamp. `a`
//! still moves with contiguity, so both it and the choice between the two ways of setting
//! the curve are what is left to check.

use rosella::embedding::umap::{Curve, EmbedOverrides, curve_params, n_x};

const LENGTHS: [usize; 20] = [
    1500, 1800, 2400, 3100, 4200, 5600, 7000, 9500, 12000, 15000, 21000, 34000, 58000, 90000,
    150000, 260000, 480000, 750000, 1200000, 2100000,
];

#[test]
fn n_x_matches_flight() {
    for (percent, expected) in [
        (10.0, 260000),
        (25.0, 750000),
        (50.0, 1200000),
        (75.0, 2100000),
        (90.0, 2100000),
    ] {
        assert_eq!(n_x(&LENGTHS, percent), expected, "at {}%", percent);
    }
    assert_eq!(n_x(&[], 50.0), 0);
}

/// `a` has to keep moving with contiguity, because the reason `b` was dropped was that it
/// could not. A fragmented assembly floors it; a contiguous one does not.
#[test]
fn a_follows_contiguity_and_b_does_not() {
    let contiguous = curve_params(&LENGTHS);
    let fragmented = curve_params(&vec![1500; 50]);

    assert!(
        contiguous.a > fragmented.a,
        "{} did not beat {}",
        contiguous.a,
        fragmented.a
    );
    assert!(
        (fragmented.a - 1.4).abs() < 1e-6,
        "floor was {}",
        fragmented.a
    );
    assert_eq!(contiguous.b, fragmented.b);
}

#[test]
fn the_curve_is_pinned_until_min_dist_or_spread_asks_for_a_fit() {
    let pinned = |overrides: &EmbedOverrides| {
        matches!(Curve::from_overrides(&LENGTHS, overrides), Curve::Pinned(_))
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
    let Curve::Fit { min_dist, spread } = Curve::from_overrides(
        &LENGTHS,
        &EmbedOverrides {
            spread: Some(2.0),
            ..Default::default()
        },
    ) else {
        panic!("spread did not ask for a fit");
    };
    assert_eq!((min_dist, spread), (0.0, 2.0));
}
