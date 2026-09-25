use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::coverage::scatter::{Scatter, fit, fit_neighbours};

const BINS: usize = 150;
const MEMBERS: usize = 40;
const PASSENGER_SHARE: f64 = 0.1;

fn normal(rng: &mut StdRng) -> f64 {
    let (u, v) = (rng.random::<f64>().max(1e-12), rng.random::<f64>());
    (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos()
}

fn simulate(truth: Scatter, rng: &mut StdRng) -> (Array2<f64>, Vec<usize>, Vec<Vec<usize>>) {
    let centres = (0..BINS)
        .map(|_| 10f64.powf(rng.random_range(-0.5..2.0)))
        .collect::<Vec<_>>();
    let mut table = Array2::zeros((BINS * MEMBERS, 2));
    let mut lengths = Vec::new();
    let mut bins = vec![Vec::new(); BINS];
    for contig in 0..BINS * MEMBERS {
        let bin = contig / MEMBERS;
        let length = rng.random_range(2_000..50_000);
        let home = match rng.random::<f64>() < PASSENGER_SHARE {
            true => rng.random_range(0..BINS),
            false => bin,
        };
        let centre = centres[home];
        let spread =
            (truth.sampling * centre / length as f64 + truth.bias * centre * centre).sqrt();
        table[[contig, 0]] = (centre + spread * normal(rng)).max(0.0);
        lengths.push(length);
        bins[bin].push(contig);
    }
    (table, lengths, bins)
}

#[test]
fn recovers_the_scatter_through_passengers() {
    let mut rng = StdRng::seed_from_u64(11);
    for truth in [
        Scatter {
            sampling: 150.0,
            bias: 0.0,
        },
        Scatter {
            sampling: 150.0,
            bias: 0.01,
        },
        Scatter {
            sampling: 3_000.0,
            bias: 0.04,
        },
    ] {
        let (table, lengths, bins) = simulate(truth, &mut rng);
        let model = fit(&table, &lengths, &bins)[0].expect("enough contigs to fit");
        for (depth, length) in [(0.5, 3_000), (5.0, 10_000), (50.0, 30_000)] {
            let ratio = model.variance(depth, length) / truth.variance(depth, length);
            assert!(
                (0.5..2.0).contains(&ratio),
                "{truth:?} fitted as {model:?}, variance off by {ratio:.2} at depth {depth}"
            );
        }
    }
}

#[test]
fn recovers_the_scatter_from_neighbours_that_are_half_strangers() {
    const K: usize = 10;
    let mut rng = StdRng::seed_from_u64(5);
    let truth = Scatter {
        sampling: 300.0,
        bias: 0.02,
    };
    let (table, lengths, _) = simulate(truth, &mut rng);
    let pool = (0..BINS * MEMBERS).collect::<Vec<_>>();
    let neighbours = pool
        .iter()
        .flat_map(|contig| {
            let genome = contig / MEMBERS;
            (0..K)
                .map(|at| match at % 2 {
                    0 => genome * MEMBERS + (contig + 1 + at) % MEMBERS,
                    _ => (contig + 1 + 97 * at) % (BINS * MEMBERS),
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let model = fit_neighbours(&table, &lengths, &pool, &neighbours, K, &mut rng)[0]
        .expect("enough pairs to fit");
    for (depth, length) in [(0.5, 3_000), (5.0, 10_000), (50.0, 30_000)] {
        let ratio = model.variance(depth, length) / truth.variance(depth, length);
        assert!(
            (0.5..2.0).contains(&ratio),
            "fitted as {model:?}, variance off by {ratio:.2} at depth {depth}"
        );
    }
}
