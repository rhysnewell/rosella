use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::recover::homing::{Sample, chances};

fn draw(rng: &mut StdRng, low: f64, high: f64) -> Sample {
    Sample {
        share: rng.random_range(low..high),
        length: rng.random_range(300..1500),
    }
}

#[test]
fn the_fitted_mix_recovers_how_many_contigs_have_a_home() {
    let mut rng = StdRng::seed_from_u64(7);
    let homed = (0..2000)
        .map(|_| (draw(&mut rng, 0.55, 1.0), true))
        .collect::<Vec<_>>();
    let homeless = (0..2000)
        .map(|_| draw(&mut rng, 0.0, 0.6))
        .collect::<Vec<_>>();
    for home in [0.2, 0.7] {
        let with = (4000.0 * home) as usize;
        let real = (0..4000)
            .map(|at| match at < with {
                true => draw(&mut rng, 0.55, 1.0),
                false => draw(&mut rng, 0.0, 0.6),
            })
            .collect::<Vec<_>>();
        let fit = chances(&homed, &homeless, &real).expect("every input is populated");
        let mean = fit.chances.iter().sum::<f64>() / fit.chances.len() as f64;
        assert!(
            (mean - home).abs() < 0.05,
            "mean chance {mean} against {home}"
        );
        let clear_but_attached = real[with..]
            .iter()
            .zip(&fit.chances[with..])
            .filter(|(sample, chance)| sample.share < 0.45 && **chance > 0.5)
            .count();
        assert_eq!(clear_but_attached, 0);
    }
}
