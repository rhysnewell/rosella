use rosella::quality::edge_spread;

type Side = (f64, Vec<(usize, f64)>);

fn worth_sum(sides: &[Side], weights: &[f64]) -> f64 {
    sides
        .iter()
        .map(|(sign, points)| {
            let total = points
                .iter()
                .map(|(marker, _)| weights[*marker])
                .sum::<f64>();
            let worth = points
                .iter()
                .map(|(marker, point)| weights[*marker] * point)
                .sum::<f64>()
                / total;
            sign * worth.max(0.0).powi(2)
        })
        .sum()
}

#[test]
fn the_spread_is_the_gradient_over_marker_weights() {
    let sides = vec![
        (1.0, vec![(0, 100.0), (1, 100.0), (2, 0.0), (3, 100.0)]),
        (1.0, vec![(2, 100.0), (4, -100.0), (5, 0.0)]),
        (-1.0, vec![(0, 100.0), (1, 0.0), (2, 100.0), (3, 100.0)]),
        (-1.0, vec![(1, 100.0), (4, 100.0), (5, 0.0)]),
        (-1.0, vec![(3, -100.0), (5, 0.0)]),
    ];
    let step = 1e-6;
    let slopes = (0..6).map(|marker| {
        let mut up = vec![1.0; 6];
        let mut down = vec![1.0; 6];
        up[marker] += step;
        down[marker] -= step;
        (worth_sum(&sides, &up) - worth_sum(&sides, &down)) / (2.0 * step)
    });
    let expected = slopes.map(|slope| slope * slope).sum::<f64>().sqrt();
    let spread = edge_spread(sides);
    assert!(
        (spread - expected).abs() < 1e-5 * expected,
        "{spread} {expected}"
    );
}

#[test]
fn two_clean_halves_joined_have_no_spread() {
    let whole = (0..8).map(|marker| (marker, 100.0)).collect::<Vec<_>>();
    let half = |range: std::ops::Range<usize>| {
        (0..8)
            .map(|marker| (marker, if range.contains(&marker) { 100.0 } else { 0.0 }))
            .collect::<Vec<_>>()
    };
    let spread = edge_spread([(1.0, whole), (-1.0, half(0..4)), (-1.0, half(4..8))]);
    assert!(spread.abs() < 1e-9, "{spread}");
}
