use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::embedding::knn::Metric;
use rosella::embedding::metrics::{DistanceSettings, MIN_VAR, prepared::PreparedAggregate};

const ROWS: usize = 41;

fn prepared(samples: usize, aggregate_weight: Option<f64>, calibrate: bool) -> PreparedAggregate {
    let mut rng = StdRng::seed_from_u64(7);
    let coverage = Array2::from_shape_fn((ROWS, 2 * samples), |(_, column)| match column % 2 {
        0 => rng.random_range(0.0..30.0),
        _ => rng.random_range(0.0..5.0),
    });
    let tnf = Array2::from_shape_fn((ROWS, 136), |(row, _)| match row {
        0 => 0.0,
        _ => rng.random_range(-2.0..2.0),
    });
    let indices = (0..ROWS).collect::<Vec<_>>();
    let lengths = (0..ROWS).map(|row| 1000 + 97 * row).collect::<Vec<_>>();
    let settings = DistanceSettings {
        presence_fraction: 0.1,
        aggregate_weight,
        calibrate,
    };
    PreparedAggregate::new(
        &coverage,
        &tnf,
        &indices,
        &[MIN_VAR; ROWS],
        &lengths,
        settings,
    )
}

// A lane that summed in another order would move neighbours by a bit and fail nothing else.
#[test]
fn a_batch_measures_every_pair_as_it_would_alone() {
    for (aggregate_weight, calibrate) in [(None, false), (Some(0.0), false), (Some(0.6), true)] {
        let prepared = prepared(2, aggregate_weight, calibrate);
        for a in [0, 5, 40] {
            let others = (0..ROWS as u32)
                .filter(|b| *b as usize != a)
                .collect::<Vec<_>>();
            let mut batch = Vec::new();
            (&prepared).within_many(a, &others, |_| f64::INFINITY, &mut batch);
            for (b, measured) in others.iter().zip(batch) {
                assert_eq!(
                    measured.to_bits(),
                    prepared.distance(a, *b as usize).to_bits()
                );
            }
        }
    }
}

// A bound that cut a pair nearer than itself would drop a true neighbour without a trace.
#[test]
fn a_bounded_pair_is_exact_or_past_its_bound() {
    for aggregate_weight in [None, Some(0.6)] {
        let prepared = prepared(6, aggregate_weight, false);
        for a in 0..ROWS {
            for b in (0..ROWS).filter(|b| *b != a) {
                let exact = prepared.distance(a, b);
                let bounds = [0.05, 0.2, 0.4, 0.6, exact, exact.next_down()];
                for bound in bounds {
                    let measured = (&prepared).within(a, b, bound);
                    match exact <= bound {
                        true => assert_eq!(measured.to_bits(), exact.to_bits()),
                        false => assert!(measured > bound, "{a} {b} {bound} {measured}"),
                    }
                }
            }
        }
    }
}
