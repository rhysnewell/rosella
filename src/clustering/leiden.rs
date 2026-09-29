use crate::clustering::graph_partition::{Incident, compact, edge_weight_total, visit_order};
use crate::embedding::{Graph, row_of};

const MAX_LEVELS: usize = 20;

// One flat run of edges, so a gather walks memory in order instead of chasing a row per node.
pub(crate) struct Level {
    offsets: Vec<usize>,
    targets: Vec<u32>,
    weights: Vec<f64>,
    pub(crate) size: Vec<f64>,
}

impl Level {
    pub(crate) fn from_graph(graph: &Graph) -> Self {
        let mut level = Self::empty(vec![1.0; graph.rows()]);
        for row in 0..graph.rows() {
            let (targets, weights) = row_of(graph, row);
            level.close(
                targets
                    .iter()
                    .zip(weights)
                    .filter(|(target, _)| **target as usize != row)
                    .map(|(target, weight)| (*target as usize, *weight as f64)),
            );
        }
        level
    }

    fn empty(size: Vec<f64>) -> Self {
        Self {
            offsets: vec![0],
            targets: Vec::new(),
            weights: Vec::new(),
            size,
        }
    }

    fn close(&mut self, edges: impl Iterator<Item = (usize, f64)>) {
        for (target, weight) in edges {
            self.targets.push(target as u32);
            self.weights.push(weight);
        }
        self.offsets.push(self.targets.len());
    }

    pub(crate) fn neighbours(&self, node: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let span = self.offsets[node]..self.offsets[node + 1];
        self.targets[span.clone()]
            .iter()
            .zip(&self.weights[span])
            .map(|(target, weight)| (*target as usize, *weight))
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
            self.neighbours(node)
                .map(|(target, weight)| (of[target], weight)),
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
            for (target, _) in level.neighbours(node) {
                if of[target] != best && !queued[target] {
                    queued[target] = true;
                    queue.push_back(target);
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
        let outward = level
            .neighbours(node)
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
    let mut next = Level::empty(size);
    for (id, nodes) in members.iter().enumerate() {
        incident.gather(nodes.iter().flat_map(|node| {
            level
                .neighbours(*node)
                .map(|(target, weight)| (ids[target], weight))
        }));
        next.close(incident.iter().filter(|(to, _)| *to != id));
    }

    (next, ids)
}

pub fn leiden(graph: &Graph, sizes: Option<&[f64]>, gamma: f64, seed: u64) -> Vec<i32> {
    leiden_from(&base_level(graph, sizes), gamma, seed)
}

pub(crate) fn base_level(graph: &Graph, sizes: Option<&[f64]>) -> Level {
    let level = Level::from_graph(graph);
    match sizes {
        Some(sizes) => level.with_size(sizes.to_vec()),
        None => level,
    }
}

// Every rung starts from the same first level, so the rungs share it rather than each holding
// a copy of the whole graph.
pub(crate) fn leiden_from(base: &Level, gamma: f64, seed: u64) -> Vec<i32> {
    let mut aggregated: Option<Level> = None;
    let mut membership = (0..base.len()).collect::<Vec<_>>();
    let mut start: Option<Vec<usize>> = None;
    let mut labels = vec![0i32; base.len()];

    for _ in 0..MAX_LEVELS {
        let level = aggregated.as_ref().unwrap_or(base);
        let of = local_move(level, gamma, seed, start.as_deref());
        for (node, slot) in labels.iter_mut().enumerate() {
            *slot = of[membership[node]] as i32;
        }

        let refined = refine(level, &of, gamma, seed);
        let (next, ids) = aggregate(level, &refined);
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
        aggregated = Some(next);
    }

    compact(&labels)
}

pub fn resolutions(graph: &Graph, sizes: Option<&[f64]>, steps: usize) -> Vec<f64> {
    let edges = edge_weight_total(graph);
    let total = sizes
        .map_or(graph.rows() as f64, |sizes| sizes.iter().sum::<f64>())
        .max(1.0);
    let mean_degree = 2.0 * edges / total;
    let largest = (total / crate::tuning::LADDER_COARSEST).max(2.0);
    let smallest = (total / crate::tuning::LADDER_FINEST).max(2.0);
    if steps < 2 || largest <= smallest {
        return vec![mean_degree / largest];
    }
    let ratio = (smallest / largest).powf(1.0 / (steps - 1) as f64);
    (0..steps)
        .map(|step| mean_degree / (largest * ratio.powi(step as i32)))
        .collect()
}
