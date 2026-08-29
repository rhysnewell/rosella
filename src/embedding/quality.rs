use std::collections::HashSet;

use ndarray::Array2;

use crate::embedding::{
    knn::{KnnGraph, build_knn},
    metrics::euclidean,
};

/// How much of each contig's neighbourhood in the source metric survived the layout. Label
/// free, so it is the only thing that can rank an embedding on the assembly being binned
/// rather than on a dataset with a gold standard.
pub fn neighbour_preservation(
    embedding: &Array2<f64>,
    source: &KnnGraph,
    seed: u64,
) -> Option<f64> {
    let neighbours = source.indices.ncols();
    if neighbours == 0 || embedding.nrows() < 2 || embedding.nrows() != source.n_points() {
        return None;
    }

    let _timer = crate::timing::scope("preservation");
    let rows = embedding
        .rows()
        .into_iter()
        .map(|row| row.to_vec())
        .collect::<Vec<_>>();
    let embedded = build_knn(rows.len(), neighbours, seed, |a, b| {
        euclidean(&rows[a], &rows[b])
    });

    let mut total = 0.0;
    let mut counted = 0usize;
    for row in 0..rows.len() {
        let expected = source
            .indices
            .row(row)
            .iter()
            .copied()
            .filter(|index| *index != u32::MAX)
            .collect::<HashSet<_>>();
        if expected.is_empty() {
            continue;
        }
        let kept = embedded
            .indices
            .row(row)
            .iter()
            .filter(|index| **index != u32::MAX && expected.contains(index))
            .count();
        total += kept as f64 / expected.len() as f64;
        counted += 1;
    }

    (counted > 0).then(|| total / counted as f64)
}
