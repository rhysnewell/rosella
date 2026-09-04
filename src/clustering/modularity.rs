use crate::clustering::graph_partition::node_degrees;
use crate::embedding::{Graph, intersect::row_of};

const NO_COMMUNITIES: f64 = -1.0;

pub fn modularity(graph: &Graph, labels: &[i32], gamma: f64) -> f64 {
    let nodes = graph.rows();
    if nodes == 0 || labels.len() < nodes {
        return NO_COMMUNITIES;
    }

    let highest = labels[..nodes].iter().copied().max().unwrap_or(-1);
    if highest < 0 {
        return NO_COMMUNITIES;
    }

    let degrees = node_degrees(graph);

    let two_m = degrees.iter().sum::<f64>();
    if two_m <= 0.0 {
        return NO_COMMUNITIES;
    }

    let mut internal = vec![0.0f64; highest as usize + 1];
    let mut incident = vec![0.0f64; highest as usize + 1];

    for row in 0..nodes {
        let label = labels[row];
        if label < 0 {
            continue;
        }
        incident[label as usize] += degrees[row];

        let (targets, weights) = row_of(graph, row);
        for (target, weight) in targets.iter().zip(weights) {
            let target = *target as usize;
            if target != row && labels[target] == label {
                internal[label as usize] += *weight as f64;
            }
        }
    }

    internal
        .iter()
        .zip(&incident)
        .map(|(inside, total)| inside / two_m - gamma * (total / two_m).powi(2))
        .sum()
}
