use crate::clustering::graph_partition::node_degrees;
use crate::embedding::{Graph, intersect::row_of};

const NO_COMMUNITIES: f64 = -1.0;

pub(crate) fn plogp(value: f64) -> f64 {
    if value > 0.0 {
        value * value.log2()
    } else {
        0.0
    }
}

fn codelengths(graph: &Graph, labels: &[i32]) -> Option<(f64, f64)> {
    let nodes = graph.rows();
    if nodes == 0 || labels.len() < nodes {
        return None;
    }

    let highest = labels[..nodes].iter().copied().max().unwrap_or(-1);
    if highest < 0 {
        return None;
    }

    let degrees = node_degrees(graph);
    let two_m = degrees.iter().sum::<f64>();
    if two_m <= 0.0 {
        return None;
    }

    let mut exit = vec![0.0f64; highest as usize + 1];
    let mut inside = vec![0.0f64; highest as usize + 1];
    let mut visit_plogp = 0.0;

    for row in 0..nodes {
        let label = labels[row];
        if label < 0 {
            continue;
        }

        let visit = degrees[row] / two_m;
        visit_plogp += plogp(visit);
        inside[label as usize] += visit;

        let (targets, weights) = row_of(graph, row);
        for (target, weight) in targets.iter().zip(weights) {
            let target = *target as usize;
            if target != row && labels[target] != label {
                exit[label as usize] += *weight as f64;
            }
        }
    }

    let exit = exit
        .iter()
        .map(|weight| weight / two_m)
        .collect::<Vec<f64>>();
    let total_exit = exit.iter().sum::<f64>();
    let modules = exit
        .iter()
        .zip(&inside)
        .map(|(out, within)| plogp(out + within) - 2.0 * plogp(*out))
        .sum::<f64>();

    Some((plogp(total_exit) - visit_plogp + modules, -visit_plogp))
}

/// Two level map equation of Rosvall & Bergstrom (2008), in bits per step, beside the share
/// of the one module codelength it saves. Lower is better on the first and higher on the
/// second, so nothing ranks on the codelength directly.
pub fn codelength_and_saving(graph: &Graph, labels: &[i32]) -> (f64, f64) {
    match codelengths(graph, labels) {
        Some((partitioned, one_module)) if one_module > 0.0 => {
            (partitioned, 1.0 - partitioned / one_module)
        }
        Some((partitioned, _)) => (partitioned, NO_COMMUNITIES),
        None => (f64::INFINITY, NO_COMMUNITIES),
    }
}

/// The share of the one module codelength the partition saves. Ranking sorts descending
/// everywhere, and a single community lands on exactly zero, as modularity does.
pub fn codelength_saving(graph: &Graph, labels: &[i32]) -> f64 {
    codelength_and_saving(graph, labels).1
}
