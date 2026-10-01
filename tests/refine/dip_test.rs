use rosella::refine::dip::{dip, exceeds_null};

const TOLERANCE: f64 = 1e-9;

/// Reference values are Hartigan's algorithm as shipped in the diptest package.
#[test]
fn the_dip_reproduces_the_reference() {
    let uniform = [
        0.6251, 0.8972, 0.7757, 0.2252, 0.3002, 0.8736, 0.0053, 0.8212, 0.7971, 0.4679, 0.303,
        0.2784, 0.2549, 0.4451, 0.5045, 0.5535, 0.9955, 0.7927, 0.6222, 0.989, 0.2153, 0.1602,
        0.6125, 0.0439, 0.0357, 0.5149, 0.4662, 0.9172, 0.6292, 0.5141, 0.4969, 0.2475, 0.0118,
        0.1924, 0.692, 0.2006, 0.3695, 0.0037, 0.83, 0.1545,
    ];
    let two_modes = [
        -3.0559, -2.9448, -2.9681, -3.6125, -2.9619, -2.3206, -3.7736, -2.5703, -2.9403, -3.3207,
        -1.9998, -2.6189, -3.5996, -2.9627, -2.7117, -3.0944, -2.6585, -3.0333, -2.6664, -2.2807,
        2.6622, 3.1016, 2.7683, 3.0636, 2.4064, 2.7103, 2.9019, 3.4494, 3.5726, 2.3382, 2.6027,
        3.3235, 2.0038, 2.7684, 2.9514, 3.6285, 3.3447, 2.8364, 2.8157, 2.8749,
    ];
    let spike = [
        0.3047, -0.0856, -0.0607, 0.0705, -0.0242, -0.0395, -0.2228, -0.0023, -0.0887, 0.2332,
        0.1306, -0.0048, 0.1337, -0.068, 0.2104, -0.0011, 0.1167, -0.2582, 0.0693, -0.3376,
        -0.4071, -0.0609, -0.18, 0.0328, 0.449, -0.1663, -0.1248, 0.0411, 0.0986, -0.0353, 5.0,
    ];
    for (values, expected) in [
        (&uniform[..], 0.058288626609442076),
        (&two_modes[..], 0.17381581530809256),
        (&spike[..], 0.037583144296205157),
    ] {
        assert!(
            (dip(values) - expected).abs() < TOLERANCE,
            "{}",
            dip(values)
        );
    }
}

#[test]
fn repeated_values_read_as_one_weighted_point() {
    let values = [
        0.7652, 0.9092, 0.1511, 0.9334, 0.0052, 0.753, 0.8105, 0.1368, 0.4189, 0.8153, 0.0143,
        0.6285,
    ];
    let weights = [3.0, 4.0, 5.0, 3.0, 1.0, 4.0, 2.0, 2.0, 1.0, 1.0, 5.0, 2.0];
    let replicated = values
        .iter()
        .zip(&weights)
        .flat_map(|(value, weight)| std::iter::repeat_n(*value, *weight as usize))
        .collect::<Vec<_>>();

    let repeated = dip(&replicated);
    assert!(
        (repeated - 0.1496713592193667).abs() < TOLERANCE,
        "{repeated}"
    );
}

/// The null is the same sample at uniform positions, so a clean second group beats every
/// draw and a smooth ramp beats none.
#[test]
fn a_second_group_exceeds_the_null_and_a_ramp_does_not() {
    let mut values = (0..30)
        .map(|i| 0.6 + 0.03 * i as f64 / 30.0)
        .collect::<Vec<_>>();
    values.extend((0..8).map(|i| 0.06 * i as f64 / 8.0));
    let ramp = (0..38)
        .map(|i| (i as f64 / 38.0).powi(2))
        .collect::<Vec<_>>();

    assert!(exceeds_null(&values, 400, 7));
    assert!(!exceeds_null(&ramp, 400, 7));
}
