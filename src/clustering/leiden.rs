use crate::clustering::graph_partition::{Incident, compact, edge_weight_total, visit_order};
use crate::embedding::{Graph, row_of};

const MAX_LEVELS: usize = 20;

pub const NULL_NAMES: [&str; 2] = ["cpm", "degree"];

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Null {
    #[default]
    Cpm,
    Degree,
}

impl Null {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "cpm" => Some(Self::Cpm),
            "degree" => Some(Self::Degree),
            _ => None,
        }
    }
}

pub(crate) struct Level {
    pub(crate) neighbours: Vec<Vec<(usize, f64)>>,
    pub(crate) size: Vec<f64>,
}

impl Level {
    pub(crate) fn from_graph(graph: &Graph) -> Self {
        let neighbours = (0..graph.rows())
            .map(|row| {
                let (targets, weights) = row_of(graph, row);
                targets
                    .iter()
                    .zip(weights)
                    .filter(|(target, _)| **target as usize != row)
                    .map(|(target, weight)| (*target as usize, *weight as f64))
                    .collect()
            })
            .collect::<Vec<Vec<_>>>();
        let size = vec![1.0; graph.rows()];
        Self { neighbours, size }
    }

    pub(crate) fn with_size(mut self, size: Vec<f64>) -> Self {
        self.size = size;
        self
    }

    pub(crate) fn len(&self) -> usize {
        self.size.len()
    }

    pub(crate) fn gather(&self, incident: &mut Incident, node: usize, of: &[usize]) {
        incident.gather(
            self.neighbours[node]
                .iter()
                .map(|(target, weight)| (of[*target], *weight)),
        );
    }
}

fn community_sizes(level: &Level, of: &[usize]) -> Vec<f64> {
    let mut sizes = vec![0.0; level.len()];
    for (node, community) in of.iter().enumerate() {
        sizes[*community] += level.size[node];
    }
    sizes
}

fn local_move(level: &Level, gamma: f64, seed: u64, start: Option<&[usize]>) -> Vec<usize> {
    let mut of = start.map_or_else(|| (0..level.len()).collect::<Vec<_>>(), <[usize]>::to_vec);
    let mut sizes = community_sizes(level, &of);
    let order = visit_order(level.len(), seed);
    let mut queued = vec![true; level.len()];
    let mut queue = order
        .iter()
        .copied()
        .collect::<std::collections::VecDeque<_>>();
    let mut incident = Incident::new(level.len());

    while let Some(node) = queue.pop_front() {
        queued[node] = false;
        let current = of[node];
        level.gather(&mut incident, node, &of);
        sizes[current] -= level.size[node];

        let mut best = current;
        let mut best_gain = incident.get(current) - gamma * level.size[node] * sizes[current];
        for (community, weight) in incident.iter() {
            let gain = weight - gamma * level.size[node] * sizes[community];
            if gain > best_gain || (gain == best_gain && community < best) {
                best = community;
                best_gain = gain;
            }
        }

        sizes[best] += level.size[node];
        if best != current {
            of[node] = best;
            for (target, _) in &level.neighbours[node] {
                if of[*target] != best && !queued[*target] {
                    queued[*target] = true;
                    queue.push_back(*target);
                }
            }
        }
    }
    of
}

/// Ties break on the lower community, as in `local_move`. Without it the winner follows
/// whatever order the incident weights happen to be visited in.
fn refine(level: &Level, of: &[usize], gamma: f64, seed: u64) -> Vec<usize> {
    let mut refined = (0..level.len()).collect::<Vec<_>>();
    let mut sizes = level.size.clone();
    let outer = community_sizes(level, of);
    let mut incident = Incident::new(level.len());

    for node in visit_order(level.len(), seed.wrapping_add(1)) {
        if sizes[refined[node]] != level.size[node] {
            continue;
        }
        let community = of[node];
        let outward = level.neighbours[node]
            .iter()
            .filter(|(target, _)| of[*target] == community)
            .map(|(_, weight)| weight)
            .sum::<f64>();
        if outward < gamma * level.size[node] * (outer[community] - level.size[node]) {
            continue;
        }

        let mut best = refined[node];
        let mut best_gain = 0.0;
        level.gather(&mut incident, node, &refined);
        for (candidate, weight) in incident.iter() {
            if candidate == refined[node] || of[candidate] != community {
                continue;
            }
            let gain = weight - gamma * level.size[node] * sizes[candidate];
            if gain > best_gain || (gain == best_gain && candidate < best) {
                best = candidate;
                best_gain = gain;
            }
        }

        if best != refined[node] {
            sizes[refined[node]] -= level.size[node];
            sizes[best] += level.size[node];
            refined[node] = best;
        }
    }
    refined
}

pub(crate) fn aggregate(level: &Level, refined: &[usize]) -> (Level, Vec<usize>) {
    let ids = compact(&refined.iter().map(|c| *c as i32).collect::<Vec<_>>())
        .into_iter()
        .map(|c| c as usize)
        .collect::<Vec<_>>();
    let count = ids.iter().copied().max().map_or(0, |m| m + 1);

    let mut size = vec![0.0; count];
    for (node, id) in ids.iter().enumerate() {
        size[*id] += level.size[node];
    }

    let mut members = vec![Vec::new(); count];
    for (node, id) in ids.iter().enumerate() {
        members[*id].push(node);
    }

    let mut incident = Incident::new(count);
    let neighbours = members
        .iter()
        .enumerate()
        .map(|(id, nodes)| {
            incident.gather(nodes.iter().flat_map(|node| {
                level.neighbours[*node]
                    .iter()
                    .map(|(target, weight)| (ids[*target], *weight))
            }));
            incident.iter().filter(|(to, _)| *to != id).collect()
        })
        .collect();

    (Level { neighbours, size }, ids)
}

pub fn leiden(graph: &Graph, sizes: Option<&[f64]>, gamma: f64, seed: u64) -> Vec<i32> {
    let mut level = Level::from_graph(graph);
    if let Some(sizes) = sizes {
        level = level.with_size(sizes.to_vec());
    }
    let mut membership = (0..graph.rows()).collect::<Vec<_>>();
    let mut start: Option<Vec<usize>> = None;
    let mut labels = vec![0i32; graph.rows()];

    for _ in 0..MAX_LEVELS {
        let of = local_move(&level, gamma, seed, start.as_deref());
        for (node, slot) in labels.iter_mut().enumerate() {
            *slot = of[membership[node]] as i32;
        }

        let refined = refine(&level, &of, gamma, seed);
        let (next, ids) = aggregate(&level, &refined);
        if next.len() == level.len() {
            break;
        }

        let mut inherited = vec![0usize; next.len()];
        for (node, id) in ids.iter().enumerate() {
            inherited[*id] = of[node];
        }
        start = Some(
            compact(&inherited.iter().map(|c| *c as i32).collect::<Vec<_>>())
                .into_iter()
                .map(|c| c as usize)
                .collect(),
        );

        for slot in &mut membership {
            *slot = ids[*slot];
        }
        level = next;
    }

    compact(&labels)
}

/// A rung aims at a community mass. `band` names that mass range directly, in the same bases
/// the sizes are in; without one the ladder divides the assembly's own mass, which puts the
/// whole ladder somewhere else on every assembly.
pub fn resolutions(
    graph: &Graph,
    sizes: Option<&[f64]>,
    steps: usize,
    band: Option<(f64, f64)>,
) -> Vec<f64> {
    let edges = edge_weight_total(graph);
    let total = sizes
        .map_or(graph.rows() as f64, |sizes| sizes.iter().sum::<f64>())
        .max(1.0);
    let mean_degree = 2.0 * edges / total;
    let (largest, smallest) = match band {
        Some((floor, ceiling)) => (ceiling.min(total).max(2.0), floor.min(total).max(2.0)),
        None => (
            (total / crate::tuning::LADDER_COARSEST).max(2.0),
            (total / crate::tuning::LADDER_FINEST).max(2.0),
        ),
    };
    if steps < 2 || largest <= smallest {
        return vec![mean_degree / largest];
    }
    let ratio = (smallest / largest).powf(1.0 / (steps - 1) as f64);
    (0..steps)
        .map(|step| mean_degree / (largest * ratio.powi(step as i32)))
        .collect()
}
