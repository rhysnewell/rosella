//! The k-NN graph is the one stage that has to be reproducible for the whole pipeline to
//! be, so it is checked for determinism as well as accuracy.

use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::embedding::knn::{KnnGraph, MAX_CANDIDATES, build_knn_with};
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

#[test]
fn repeated_builds_agree() {
    let rows = sample_rows(400, 8, 11);

    let first = build_knn_with(rows.len(), 15, MAX_CANDIDATES, 42, |i, j| euclidean(&rows[i], &rows[j]));
    let second = build_knn_with(rows.len(), 15, MAX_CANDIDATES, 42, |i, j| euclidean(&rows[i], &rows[j]));

    assert_eq!(first.indices, second.indices);
    assert_eq!(first.dists, second.dists);
}

#[test]
fn a_different_seed_still_finds_the_same_neighbours() {
    let rows = sample_rows(400, 8, 11);
    let exact = exact_knn(&rows, 15);

    let first = build_knn_with(rows.len(), 15, MAX_CANDIDATES, 1, |i, j| euclidean(&rows[i], &rows[j]));
    let second = build_knn_with(rows.len(), 15, MAX_CANDIDATES, 99999, |i, j| euclidean(&rows[i], &rows[j]));

    for row in 0..rows.len() {
        let exact_row: Vec<u32> = exact.indices.row(row).to_vec();
        assert!(recall(&first.indices.row(row).to_vec(), &exact_row) > 0.8);
        assert!(recall(&second.indices.row(row).to_vec(), &exact_row) > 0.8);
    }
}

#[test]
fn descent_recovers_the_exact_neighbours() {
    let rows = sample_rows(500, 6, 3);

    let approximate = build_knn_with(rows.len(), 10, MAX_CANDIDATES, 42, |i, j| euclidean(&rows[i], &rows[j]));
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
    let graph = build_knn_with(rows.len(), 12, MAX_CANDIDATES, 42, |i, j| euclidean(&rows[i], &rows[j]));

    for row in 0..rows.len() {
        let indices = graph.indices.row(row);
        let dists = graph.dists.row(row);
        assert!(!indices.iter().any(|index| *index as usize == row));
        assert!(dists.windows(2).into_iter().all(|pair| pair[0] <= pair[1]));
    }
}
