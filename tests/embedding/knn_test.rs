//! The k-NN graph is the one stage that has to be reproducible for the whole pipeline to
//! be, so it is checked for determinism as well as accuracy.

use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::embedding::knn::{KnnGraph, build_knn_from, build_knn_with, candidates, nearest_in};
use rosella::embedding::metrics::euclidean;

fn exact_knn(rows: &[Vec<f64>], k: usize) -> KnnGraph {
    let mut indices = Array2::from_elem((rows.len(), k), u32::MAX);
    let mut dists = Array2::from_elem((rows.len(), k), f32::INFINITY);
    for i in 0..rows.len() {
        let mut ranked = (0..rows.len())
            .filter(|j| *j != i)
            .map(|j| (euclidean(&rows[i], &rows[j]), j as u32))
            .collect::<Vec<_>>();
        ranked.sort_by(|a, b| a.partial_cmp(b).expect("distances are finite"));
        for (slot, (distance, neighbour)) in ranked.into_iter().take(k).enumerate() {
            indices[[i, slot]] = neighbour;
            dists[[i, slot]] = distance as f32;
        }
    }
    KnnGraph { indices, dists }
}

fn sample_rows(n_points: usize, n_features: usize, seed: u64) -> Vec<Vec<f64>> {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..n_points)
        .map(|_| {
            (0..n_features)
                .map(|_| rng.random_range(-5.0..5.0))
                .collect()
        })
        .collect()
}

fn recall(approximate: &[u32], exact: &[u32]) -> f64 {
    let found = approximate
        .iter()
        .filter(|index| exact.contains(index))
        .count();
    found as f64 / exact.len() as f64
}

/// These exercise the descent at a k far below the shipped 100, where half of k is a handful
/// of candidates. They hand it a wide cap so the subject is the descent and not the rate.
fn wide(k: usize) -> usize {
    2 * k
}

fn build_on(threads: usize, rows: &[Vec<f64>], k: usize) -> KnnGraph {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("a thread pool")
        .install(|| {
            build_knn_with(rows.len(), k, wide(k), 42, |i, j| {
                euclidean(&rows[i], &rows[j])
            })
        })
}

/// Two builds at one thread count agree for reasons unrelated to the descent, so the pool
/// width is what varies, over an instance large enough to take several rounds.
#[test]
fn builds_agree_whatever_the_thread_count() {
    let rows = sample_rows(6000, 12, 11);

    let narrow = build_on(1, &rows, 25);
    let wide = build_on(8, &rows, 25);

    assert_eq!(narrow.indices, wide.indices);
    assert_eq!(narrow.dists, wide.dists);
}

/// The shipped graph is k=100, where the rate gives 50 candidates. A thin cap costs recall
/// against the exact neighbours, so the number the binner actually runs at is the one to pin.
#[test]
fn the_shipped_candidate_rate_recovers_the_neighbours() {
    let rows = sample_rows(600, 8, 3);
    let k = 100;
    let exact = exact_knn(&rows, k);
    let built = build_knn_with(rows.len(), k, candidates(k), 42, |i, j| {
        euclidean(&rows[i], &rows[j])
    });

    let mean = (0..rows.len())
        .map(|row| {
            recall(
                &built.indices.row(row).to_vec(),
                &exact.indices.row(row).to_vec(),
            )
        })
        .sum::<f64>()
        / rows.len() as f64;
    assert!(mean > 0.95, "mean recall at the shipped rate is {mean}");
}

#[test]
fn a_different_seed_still_finds_the_same_neighbours() {
    let rows = sample_rows(400, 8, 11);
    let exact = exact_knn(&rows, 15);

    let first = build_knn_with(rows.len(), 15, wide(15), 1, |i, j| {
        euclidean(&rows[i], &rows[j])
    });
    let second = build_knn_with(rows.len(), 15, wide(15), 99999, |i, j| {
        euclidean(&rows[i], &rows[j])
    });

    for row in 0..rows.len() {
        let exact_row: Vec<u32> = exact.indices.row(row).to_vec();
        assert!(recall(&first.indices.row(row).to_vec(), &exact_row) > 0.8);
        assert!(recall(&second.indices.row(row).to_vec(), &exact_row) > 0.8);
    }
}

#[test]
fn descent_recovers_the_exact_neighbours() {
    let rows = sample_rows(500, 6, 3);

    let approximate = build_knn_with(rows.len(), 10, wide(10), 42, |i, j| {
        euclidean(&rows[i], &rows[j])
    });
    let exact = exact_knn(&rows, 10);

    let mean_recall = (0..rows.len())
        .map(|row| {
            recall(
                &approximate.indices.row(row).to_vec(),
                &exact.indices.row(row).to_vec(),
            )
        })
        .sum::<f64>()
        / rows.len() as f64;

    assert!(mean_recall > 0.95, "mean recall was {}", mean_recall);
}

#[test]
fn neighbours_are_sorted_and_exclude_self() {
    let rows = sample_rows(200, 4, 7);
    let graph = build_knn_with(rows.len(), 12, wide(12), 42, |i, j| {
        euclidean(&rows[i], &rows[j])
    });

    for row in 0..rows.len() {
        let indices = graph.indices.row(row);
        let dists = graph.dists.row(row);
        assert!(!indices.iter().any(|index| *index as usize == row));
        assert!(dists.windows(2).into_iter().all(|pair| pair[0] <= pair[1]));
    }
}

fn exact_in(base: &[Vec<f64>], queries: &[Vec<f64>], k: usize) -> Vec<Vec<u32>> {
    queries
        .iter()
        .map(|query| {
            let mut ranked = (0..base.len())
                .map(|at| (euclidean(query, &base[at]), at as u32))
                .collect::<Vec<_>>();
            ranked.sort_by(|a, b| a.partial_cmp(b).expect("distances are finite"));
            ranked.into_iter().take(k).map(|(_, at)| at).collect()
        })
        .collect()
}

/// Queries drawn from the base's own spread have to be walked to from random starts, at the rate
/// the binner runs, and the answer must not hang on how many threads searched.
#[test]
fn a_query_walk_finds_its_nearest_base_points_on_any_pool() {
    let base = sample_rows(1500, 8, 5);
    let queries = sample_rows(200, 8, 6);
    let k = 30;
    let graph = build_knn_with(base.len(), k, candidates(k), 42, |i, j| {
        euclidean(&base[i], &base[j])
    });
    let search = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("a thread pool")
            .install(|| {
                nearest_in(&graph, queries.len(), k, candidates(k), 42, |query, at| {
                    euclidean(&queries[query], &base[at])
                })
            })
    };
    let narrow = search(1);
    let wide = search(8);
    assert_eq!(narrow.indices, wide.indices);

    let exact = exact_in(&base, &queries, k);
    let mean = (0..queries.len())
        .map(|query| recall(&narrow.indices.row(query).to_vec(), &exact[query]))
        .sum::<f64>()
        / queries.len() as f64;
    assert!(mean > 0.95, "mean query recall is {mean}");
}

/// The pass grows a settled graph by rows whose nearest old rows are known. Starting from those
/// lists has to reach the neighbours a descent from random reaches.
#[test]
fn a_seeded_descent_reaches_the_exact_neighbours_of_the_grown_set() {
    let rows = sample_rows(900, 8, 21);
    let old = 600;
    let k = 30;
    let exact = exact_knn(&rows, k);
    let within = exact_in(&rows[..old], &rows, k);
    let start = Array2::from_shape_fn((rows.len(), k), |(row, slot)| {
        match within[row]
            .iter()
            .filter(|at| **at as usize != row)
            .nth(slot)
        {
            Some(at) => *at,
            None => u32::MAX,
        }
    });
    let built = build_knn_from(&start, k, candidates(k), 42, |i, j| {
        euclidean(&rows[i], &rows[j])
    });
    let mean = (0..rows.len())
        .map(|row| {
            recall(
                &built.indices.row(row).to_vec(),
                &exact.indices.row(row).to_vec(),
            )
        })
        .sum::<f64>()
        / rows.len() as f64;
    assert!(mean > 0.95, "mean recall from the seeded lists is {mean}");
}
