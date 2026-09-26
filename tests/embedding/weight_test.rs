use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::embedding::weight::{Contigs, centre, recall};

const GENOMES: usize = 30;
const PER_GENOME: usize = 8;
const SAMPLES: usize = 3;
const WIDTH: usize = 32;

struct Spread {
    composition_between: f64,
    half_noise: f64,
    depth_between: f64,
    depth_within: f64,
}

fn normal(rng: &mut StdRng) -> f64 {
    let (u, v) = (rng.random::<f64>().max(1e-12), rng.random::<f64>());
    (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos()
}

fn weight_for(spread: Spread) -> f64 {
    let mut rng = StdRng::seed_from_u64(7);
    let mut coverage = Vec::new();
    let mut whole = Vec::new();
    let mut first = Vec::new();
    let mut second = Vec::new();
    for _ in 0..GENOMES {
        let centre = (0..WIDTH)
            .map(|_| spread.composition_between * normal(&mut rng))
            .collect::<Vec<_>>();
        let depths = (0..SAMPLES)
            .map(|_| (3.0 + spread.depth_between * normal(&mut rng)).exp())
            .collect::<Vec<_>>();
        for _ in 0..PER_GENOME {
            let contig = centre
                .iter()
                .map(|value| value + 0.02 * normal(&mut rng))
                .collect::<Vec<_>>();
            let mut half = || {
                contig
                    .iter()
                    .map(|value| value + spread.half_noise * normal(&mut rng))
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
        first,
        second,
    };
    centre(&recall(&contigs, 0.01, 42).expect("enough contigs to judge"))
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
            },
            false,
        ),
    ];
    for (name, spread, coverage_wins) in cases {
        let weight = weight_for(spread);
        assert_eq!(weight > 0.5, coverage_wins, "{name}: weight {weight:.3}");
    }
}
