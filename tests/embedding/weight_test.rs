use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::embedding::metrics::DistanceSettings;
use rosella::embedding::weight::{Contigs, Derived, derive};

const GENOMES: usize = 30;
const PER_GENOME: usize = 8;
const SAMPLES: usize = 3;
const WIDTH: usize = 32;

struct Spread {
    composition_between: f64,
    half_noise: f64,
    depth_between: f64,
    depth_within: f64,
    lengths_vary: bool,
}

fn normal(rng: &mut StdRng) -> f64 {
    let (u, v) = (rng.random::<f64>().max(1e-12), rng.random::<f64>());
    (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos()
}

fn derived_for(spread: Spread) -> Derived {
    let mut rng = StdRng::seed_from_u64(7);
    let mut coverage = Vec::new();
    let mut whole = Vec::new();
    let mut first = Vec::new();
    let mut second = Vec::new();
    let mut lengths = Vec::new();
    for _ in 0..GENOMES {
        let centre = (0..WIDTH)
            .map(|_| spread.composition_between * normal(&mut rng))
            .collect::<Vec<_>>();
        let depths = (0..SAMPLES)
            .map(|_| (3.0 + spread.depth_between * normal(&mut rng)).exp())
            .collect::<Vec<_>>();
        for _ in 0..PER_GENOME {
            let length = match spread.lengths_vary {
                true => 10f64.powf(rng.random_range(3.5..5.0)) as usize,
                false => 10_000,
            };
            let noise = spread.half_noise * (10_000.0 / length as f64).sqrt();
            lengths.push(length);
            let contig = centre
                .iter()
                .map(|value| value + 0.02 * normal(&mut rng))
                .collect::<Vec<_>>();
            let mut half = || {
                contig
                    .iter()
                    .map(|value| value + noise * normal(&mut rng))
                    .collect::<Vec<_>>()
            };
            first.push(half());
            second.push(half());
            whole.push(contig);
            coverage.push(
                depths
                    .iter()
                    .flat_map(|depth| {
                        let mean = depth * (spread.depth_within * normal(&mut rng)).exp();
                        [mean, mean]
                    })
                    .collect::<Vec<_>>(),
            );
        }
    }
    let contigs = Contigs {
        coverage: coverage.iter().map(Vec::as_slice).collect(),
        whole: whole.iter().map(Vec::as_slice).collect(),
        lengths,
        first,
        second,
    };
    derive(
        &contigs,
        DistanceSettings {
            presence_fraction: 0.01,
            ..DistanceSettings::default()
        },
        42,
    )
    .expect("enough contigs to judge")
}

#[test]
fn the_weight_follows_the_view_that_separates_genomes() {
    let cases = [
        (
            "coverage separates",
            Spread {
                composition_between: 0.02,
                half_noise: 0.3,
                depth_between: 1.0,
                depth_within: 0.05,
                lengths_vary: false,
            },
            true,
        ),
        (
            "composition separates",
            Spread {
                composition_between: 1.0,
                half_noise: 0.1,
                depth_between: 0.05,
                depth_within: 0.5,
                lengths_vary: false,
            },
            false,
        ),
    ];
    for (name, spread, coverage_wins) in cases {
        let weight = derived_for(spread).weight;
        assert_eq!(weight > 0.5, coverage_wins, "{name}: weight {weight:.3}");
    }
}

#[test]
fn short_contigs_lean_on_coverage_when_their_composition_is_noisier() {
    let line = derived_for(Spread {
        composition_between: 0.3,
        half_noise: 0.3,
        depth_between: 0.3,
        depth_within: 0.1,
        lengths_vary: true,
    })
    .line;
    assert!(line.slope < 0.0, "{line:?}");
}
