use std::collections::HashSet;

use ndarray::Array2;

// A paralog's copies sit in the bin's own core, where a second genome's copy sits at its edge.
pub fn outside(contigs: &[usize], lengths: &[usize], points: &Array2<f64>) -> HashSet<usize> {
    let placed = contigs
        .iter()
        .copied()
        .filter(|contig| *contig < points.nrows())
        .collect::<Vec<_>>();
    let weight = |contig: usize| lengths.get(contig).copied().unwrap_or(0) as f64;
    let total = placed.iter().map(|contig| weight(*contig)).sum::<f64>();
    let mut summed = vec![0.0; points.ncols()];
    for contig in &placed {
        for (sum, value) in summed.iter_mut().zip(points.row(*contig)) {
            *sum += weight(*contig) * value;
        }
    }
    let distances = placed
        .iter()
        .map(|contig| {
            let rest = total - weight(*contig);
            if rest <= 0.0 {
                return 0.0;
            }
            points
                .row(*contig)
                .iter()
                .zip(&summed)
                .map(|(value, sum)| (value - (sum - weight(*contig) * value) / rest).powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .collect::<Vec<_>>();
    let median = median(&distances);
    contigs
        .iter()
        .copied()
        .filter(|contig| *contig >= points.nrows())
        .chain(
            placed
                .iter()
                .zip(&distances)
                .filter(|(_, distance)| **distance >= median)
                .map(|(contig, _)| *contig),
        )
        .collect()
}

fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    match sorted.len() {
        0 => 0.0,
        n if n % 2 == 1 => sorted[n / 2],
        n => (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0,
    }
}
