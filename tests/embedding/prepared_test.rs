use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::embedding::knn::Metric;
use rosella::embedding::metrics::{DistanceSettings, MIN_VAR, prepared::PreparedAggregate};

// A lane that summed in another order would move neighbours by a bit and fail nothing else.
#[test]
fn a_batch_measures_every_pair_as_it_would_alone() {
    let rows = 41;
    let mut rng = StdRng::seed_from_u64(7);
    let coverage = Array2::from_shape_fn((rows, 4), |(_, column)| match column % 2 {
        0 => rng.random_range(0.0..30.0),
        _ => rng.random_range(0.0..5.0),
    });
    let tnf = Array2::from_shape_fn((rows, 136), |(row, _)| match row {
        0 => 0.0,
        _ => rng.random_range(-2.0..2.0),
    });
    let indices = (0..rows).collect::<Vec<_>>();
    let lengths = (0..rows).map(|row| 1000 + 97 * row).collect::<Vec<_>>();
    for (aggregate_weight, calibrate) in [(None, false), (Some(0.0), false), (Some(0.6), true)] {
        let settings = DistanceSettings {
            presence_fraction: 0.1,
            aggregate_weight,
            calibrate,
        };
        let prepared = PreparedAggregate::new(
            &coverage,
            &tnf,
            &indices,
            &[MIN_VAR; 41],
            &lengths,
            settings,
        );
        for a in [0, 5, 40] {
            let others = (0..rows as u32)
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
