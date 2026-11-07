use rosella::embedding::metrics::erfc::erfc;

#[test]
fn erfc_stays_within_machine_precision_of_libm_across_the_range() {
    let mut z = -8.0;
    while z < 30.0 {
        let expected = libm::erfc(z);
        assert!(
            (erfc(z) - expected).abs() <= 2.0 * f64::EPSILON,
            "erfc({z}) = {} against {expected}",
            erfc(z)
        );
        z += 1.0 / 4099.0;
    }
}

#[test]
fn erfc_keeps_nan_and_the_limits() {
    assert!(erfc(f64::NAN).is_nan());
    assert_eq!(erfc(f64::NEG_INFINITY), 2.0);
    assert_eq!(erfc(f64::INFINITY), 0.0);
}
