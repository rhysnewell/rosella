use std::collections::BTreeMap;

use crate::embedding::{features::ContigFeatures, metrics::AggregateMetric};
use crate::refine::bin_stats::{AGGREGATE, bin_stats, centroid};

/// `evaluate_outliers` re-embeds the noise pool alone, so a contig beside a good bin can never
/// rejoin it: that bin is not in the second layout. The bar is the bin's own contig spread.
pub fn recruit(
    features: &ContigFeatures,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    outliers: Vec<usize>,
    seed: u64,
) -> (Vec<usize>, usize) {
    if bins.is_empty() || outliers.is_empty() {
        return (outliers, 0);
    }

    let metric = AggregateMetric::new(features.n_samples() * 2, features.distance_settings());
    let mut targets = Vec::with_capacity(bins.len());
    for (id, indices) in bins.iter() {
        let Some(stats) = bin_stats(features, indices, seed) else {
            continue;
        };
        targets.push((*id, centroid(features, indices), stats.mean[AGGREGATE]));
    }
    if targets.is_empty() {
        return (outliers, 0);
    }

    let coverage_columns = features.n_samples() * 2;
    let mut recruited = 0;
    let mut left_over = Vec::new();
    for contig in outliers {
        let mut row = Vec::with_capacity(coverage_columns + features.tnf_row(contig).len());
        row.extend_from_slice(features.coverage_row(contig));
        row.extend_from_slice(features.tnf_row(contig));
        let floor = features.variance_floor(contig);

        let mut best: Option<(f64, usize)> = None;
        for (id, centre, spread) in targets.iter() {
            let distance = metric.distance(&row, &centre.row, floor, centre.floor);
            if distance <= *spread && best.is_none_or(|(closest, _)| distance < closest) {
                best = Some((distance, *id));
            }
        }

        match best {
            Some((_, id)) => {
                bins.entry(id).or_default().push(contig);
                recruited += 1;
            }
            None => left_over.push(contig),
        }
    }

    if recruited > 0 {
        for indices in bins.values_mut() {
            indices.sort_unstable();
        }
    }
    (left_over, recruited)
}
