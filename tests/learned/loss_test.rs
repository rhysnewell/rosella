use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::learned::loss::per_positive;

#[test]
fn per_positive_gradient_matches_finite_differences() {
    let mut rng = StdRng::seed_from_u64(5);
    let rows = Array2::from_shape_simple_fn((9, 4), || rng.random_range(-1.0f32..1.0));
    let labels = [0, 0, 0, 1, 1, 1, 2, 2, 2];
    let (_, analytic) = per_positive(rows.view(), &labels, 0.3);
    let step = 1e-3;
    for (row, column) in [(0, 0), (4, 2), (8, 3)] {
        let mut moved = rows.clone();
        moved[[row, column]] += step;
        let up = per_positive(moved.view(), &labels, 0.3).0;
        moved[[row, column]] -= 2.0 * step;
        let down = per_positive(moved.view(), &labels, 0.3).0;
        let numeric = (up - down) / f64::from(2.0 * step);
        let given = f64::from(analytic[[row, column]]);
        assert!(
            (numeric - given).abs() < 1e-3 + 0.02 * numeric.abs(),
            "{row},{column}: {given:.6} against {numeric:.6}"
        );
    }
}

#[test]
fn per_positive_loss_matches_a_direct_sum() {
    let mut rng = StdRng::seed_from_u64(13);
    let rows = Array2::from_shape_simple_fn((40, 8), || rng.random_range(-1.0f32..1.0));
    let labels = (0..40).map(|at| at / 5).collect::<Vec<u32>>();
    let temperature = 0.07;
    let logit = |a: usize, b: usize| {
        rows.row(a)
            .iter()
            .zip(rows.row(b))
            .map(|(x, y)| f64::from(x * y))
            .sum::<f64>()
            / f64::from(temperature)
    };
    let (mut total, mut pairs) = (0.0, 0);
    for anchor in 0..40 {
        let negatives = (0..40)
            .filter(|other| labels[*other] != labels[anchor])
            .map(|other| logit(anchor, other).exp())
            .sum::<f64>();
        for mate in (0..40).filter(|mate| *mate != anchor && labels[*mate] == labels[anchor]) {
            let own = logit(anchor, mate);
            total += (own.exp() + negatives).ln() - own;
            pairs += 1;
        }
    }
    let (given, _) = per_positive(rows.view(), &labels, temperature);
    let direct = total / f64::from(pairs);
    assert!(
        (given - direct).abs() < 1e-4 * direct,
        "{given} against {direct}"
    );
}
