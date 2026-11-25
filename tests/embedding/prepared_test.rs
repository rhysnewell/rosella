use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::embedding::knn::Metric;
use rosella::embedding::metrics::{
    DistanceSettings,
    prepared::{PreparedAggregate, QUERIES},
};

const ROWS: usize = 41;

fn prepared(samples: usize, aggregate_weight: Option<f64>) -> PreparedAggregate {
    let mut rng = StdRng::seed_from_u64(7);
    let tnf = Array2::from_shape_fn((ROWS, 136), |(row, _)| match row {
        0 => 0.0,
        _ => rng.random_range(-2.0..2.0),
    });
    prepared_on(&tnf, samples, aggregate_weight)
}

// Four tight groups, so most bounds fall between near pairs rather than past them all.
fn grouped(samples: usize, aggregate_weight: Option<f64>) -> PreparedAggregate {
    let mut rng = StdRng::seed_from_u64(11);
    let centres = Array2::from_shape_fn((4, 128), |_| rng.random_range(-1.0..1.0));
    let tnf = Array2::from_shape_fn((ROWS, 128), |(row, column)| {
        centres[[row % 4, column]] + rng.random_range(-0.2..0.2)
    });
    prepared_on(&tnf, samples, aggregate_weight)
}

fn prepared_on(
    tnf: &Array2<f64>,
    samples: usize,
    aggregate_weight: Option<f64>,
) -> PreparedAggregate {
    let mut rng = StdRng::seed_from_u64(7);
    let coverage = Array2::from_shape_fn((ROWS, 2 * samples), |(_, column)| match column % 2 {
        0 => rng.random_range(0.0..30.0),
        _ => rng.random_range(0.0..5.0),
    });
    let indices = (0..ROWS).collect::<Vec<_>>();
    let settings = DistanceSettings {
        presence_fraction: 0.1,
        aggregate_weight,
    };
    PreparedAggregate::new(&coverage, tnf, &indices, &settings)
}

// A lane that summed in another order would move neighbours by a bit and fail nothing else.
#[test]
fn a_batch_measures_every_pair_as_it_would_alone() {
    for aggregate_weight in [None, Some(0.0), Some(0.6)] {
        let prepared = prepared(2, aggregate_weight);
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
        let queries = [0, 5, 40];
        let block = prepared.block(&queries);
        let mut measured = [0.0; QUERIES];
        for b in 0..ROWS {
            prepared.within_block(&block, b, &[f64::INFINITY; QUERIES], &mut measured);
            for (a, distance) in queries.iter().zip(measured) {
                assert_eq!(distance.to_bits(), prepared.distance(*a, b).to_bits());
            }
        }
    }
}

// A bound that cut a pair nearer than itself would drop a true neighbour without a trace.
#[test]
fn a_bounded_pair_is_exact_or_past_its_bound() {
    let held = |exact: f64, bound: f64, measured: f64| match exact <= bound {
        true => measured.to_bits() == exact.to_bits(),
        false => measured > bound,
    };
    for aggregate_weight in [None, Some(0.0), Some(0.6)] {
        for prepared in [prepared(6, aggregate_weight), grouped(6, aggregate_weight)] {
            let queries = [0, 1, 2, 5, 13, 22, 39, 40];
            let block = prepared.block(&queries);
            for a in 0..ROWS {
                let others = (0..ROWS as u32)
                    .filter(|b| *b as usize != a)
                    .collect::<Vec<_>>();
                for b in others.iter().map(|b| *b as usize) {
                    let exact = prepared.distance(a, b);
                    let bounds = [0.05, 0.2, 0.4, 0.6, exact, exact.next_down()];
                    for bound in bounds {
                        let measured = (&prepared).within(a, b, bound);
                        assert!(held(exact, bound, measured), "{a} {b} {bound} {measured}");
                        let mut batch = Vec::new();
                        let walk = prepared.shifted(0);
                        walk.within_many(a, &others, |_| bound, &mut batch);
                        for (other, measured) in others.iter().zip(batch) {
                            let exact = prepared.distance(a, *other as usize);
                            assert!(held(exact, bound, measured), "{a} {other} {bound}");
                        }
                        let mut lanes = [0.0; QUERIES];
                        prepared.within_block(&block, a, &[bound; QUERIES], &mut lanes);
                        for (query, measured) in queries.iter().zip(lanes) {
                            let exact = prepared.distance(*query, a);
                            assert!(held(exact, bound, measured), "{query} {a} {bound}");
                        }
                    }
                }
            }
        }
    }
}

// The device keeps every query's exact neighbours, so it has to agree with a search of every pair.
#[test]
fn the_device_finds_the_neighbours_every_pair_gives() {
    let Some(_gpu) = rosella::learned::gemm::gpu(true) else {
        return;
    };
    let queries = (0..ROWS).collect::<Vec<_>>();
    for aggregate_weight in [Some(0.0), Some(0.6), None] {
        let prepared = grouped(6, aggregate_weight);
        // A list past 32 spans two registers of a warp, so slot 31 has to shift across lanes.
        for k in [12, ROWS - 1] {
            let device = prepared
                .searched_on_device(ROWS, &queries, k)
                .expect("an open device answers");
            let exact = rosella::embedding::knn::nearest_exact(ROWS, &queries, k, &prepared);
            for &query in &queries {
                let mut found = device.indices.row(query).to_vec();
                let mut wanted = exact.indices.row(query).to_vec();
                found.sort_unstable();
                wanted.sort_unstable();
                assert_eq!(found, wanted, "{aggregate_weight:?} {k} {query}");
                let dists = device
                    .dists
                    .row(query)
                    .into_iter()
                    .zip(exact.dists.row(query));
                for (near, far) in dists {
                    assert!((near - far).abs() < 1e-4, "{k} {query} {near} {far}");
                }
            }
        }
    }
}
