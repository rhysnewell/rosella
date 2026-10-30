use std::collections::HashMap;

use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::learned::loss::{Mates, per_positive_mated};

// Label 0 pulls label 1's first row and leaves label 2's rows out of its negatives.
fn cases() -> [Mates; 2] {
    [
        Mates::default(),
        Mates {
            pulled: HashMap::from([(0, vec![3])]),
            masked: HashMap::from([(0, vec![6, 7, 8])]),
        },
    ]
}

#[test]
fn per_positive_gradient_matches_finite_differences() {
    let mut rng = StdRng::seed_from_u64(5);
    let rows = Array2::from_shape_simple_fn((9, 4), || rng.random_range(-1.0f32..1.0));
    let labels = [0, 0, 0, 1, 1, 1, 2, 2, 2];
    for mates in cases() {
        let loss = |rows: &Array2<f32>| per_positive_mated(rows.view(), &labels, &mates, 0.3);
        let (_, analytic, _) = loss(&rows);
        let step = 1e-3;
        for (row, column) in [(0, 0), (3, 1), (4, 2), (7, 3)] {
            let mut moved = rows.clone();
            moved[[row, column]] += step;
            let up = loss(&moved).0;
            moved[[row, column]] -= 2.0 * step;
            let down = loss(&moved).0;
            let numeric = (up - down) / f64::from(2.0 * step);
            let given = f64::from(analytic[[row, column]]);
            assert!(
                (numeric - given).abs() < 1e-3 + 0.02 * numeric.abs(),
                "{row},{column}: {given:.6} against {numeric:.6}"
            );
        }
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
    let mates = Mates {
        pulled: HashMap::from([(0, vec![5, 6]), (3, vec![30])]),
        masked: HashMap::from([(0, (10..20).collect()), (7, vec![0, 1])]),
    };
    for mates in [Mates::default(), mates] {
        let list =
            |of: &HashMap<u32, Vec<u32>>, label: u32| of.get(&label).cloned().unwrap_or_default();
        let (mut total, mut pairs) = (0.0, 0);
        for anchor in 0..40 {
            let label = labels[anchor];
            let pulled = list(&mates.pulled, label);
            let masked = list(&mates.masked, label);
            let positive =
                |other: usize| labels[other] == label || pulled.contains(&(other as u32));
            let negatives = (0..40)
                .filter(|other| !positive(*other) && !masked.contains(&(*other as u32)))
                .map(|other| logit(anchor, other).exp())
                .sum::<f64>();
            for mate in (0..40).filter(|mate| *mate != anchor && positive(*mate)) {
                let own = logit(anchor, mate);
                total += (own.exp() + negatives).ln() - own;
                pairs += 1;
            }
        }
        let (given, _, _) = per_positive_mated(rows.view(), &labels, &mates, temperature);
        let direct = total / f64::from(pairs);
        assert!(
            (given - direct).abs() < 1e-4 * direct,
            "{given} against {direct}"
        );
    }
}
