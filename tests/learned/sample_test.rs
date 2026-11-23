use rand::{SeedableRng, rngs::StdRng};
use rosella::learned::sample::gamma;

#[test]
fn gamma_draws_match_the_mean_and_variance_of_their_shape() {
    let mut rng = StdRng::seed_from_u64(7);
    for shape in [1.0, 4.5, 60.0] {
        let draws = (0..200_000)
            .map(|_| gamma(&mut rng, shape))
            .collect::<Vec<_>>();
        let mean = draws.iter().sum::<f64>() / draws.len() as f64;
        let variance = draws.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / draws.len() as f64;
        assert!(
            (mean / shape - 1.0).abs() < 0.01,
            "shape {shape} mean {mean}"
        );
        assert!(
            (variance / shape - 1.0).abs() < 0.03,
            "shape {shape} variance {variance}"
        );
    }
}
