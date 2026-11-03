use std::collections::HashSet;

use ndarray::Array2;

// A paralog's copies sit in the bin's own core, where a second genome's copy sits at its edge.
// Marker carriers anchor the core, so markerless passengers cannot widen it.
pub fn outside(
    contigs: &[usize],
    lengths: &[usize],
    points: &Array2<f64>,
    carries: impl Fn(usize) -> bool,
) -> HashSet<usize> {
    let placed = contigs
        .iter()
        .copied()
        .filter(|contig| *contig < points.nrows())
        .collect::<Vec<_>>();
    let mut anchors = placed
        .iter()
        .copied()
        .filter(|contig| carries(*contig))
        .collect::<Vec<_>>();
    if anchors.len() < 2 {
        anchors.clone_from(&placed);
    }
    let weight = |contig: usize| lengths.get(contig).copied().unwrap_or(0) as f64;
    let total = anchors.iter().map(|contig| weight(*contig)).sum::<f64>();
    let mut summed = vec![0.0; points.ncols()];
    for contig in &anchors {
        for (sum, value) in summed.iter_mut().zip(points.row(*contig)) {
            *sum += weight(*contig) * value;
        }
    }
    let anchored = anchors.iter().copied().collect::<HashSet<_>>();
    let distance = |contig: usize| {
        let own = match anchored.contains(&contig) {
            true => weight(contig),
            false => 0.0,
        };
        let rest = total - own;
        if rest <= 0.0 {
            return 0.0;
        }
        points
            .row(contig)
            .iter()
            .zip(&summed)
            .map(|(value, sum)| (value - (sum - own * value) / rest).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    let median = median(
        &anchors
            .iter()
            .map(|contig| distance(*contig))
            .collect::<Vec<_>>(),
    );
    contigs
        .iter()
        .copied()
        .filter(|contig| *contig >= points.nrows() || distance(*contig) >= median)
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
