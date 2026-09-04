use rand::{Rng, SeedableRng, rngs::StdRng};

use crate::clustering::graph_partition::{Incident, Weights, compact, visit_order};
use crate::embedding::{Graph, intersect::row_of};

const MAX_LEVELS: usize = 20;

/// Far enough from the `+ 1` the visit order shuffles on that the sampler is not replaying
/// the same stream one node later.
const SAMPLE_SEED_OFFSET: u64 = 0x9E37_79B9_7F4A_7C15;

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

    /// Infomap weighs a node by its degree where CPM weighs it by one, and aggregation sums
    /// whichever it is handed.
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

/// Traag et al. (2019) draw the target with probability proportional to `exp(gain / theta)`.
/// Theta scales with the peak gain because gamma spans 256x across the ladder, so one
/// absolute value would be greedy at one end of it and uniform at the other.
fn draw(
    candidates: &mut [(usize, f64)],
    best: usize,
    peak: f64,
    theta: f64,
    rng: &mut StdRng,
) -> usize {
    if peak <= 0.0 {
        return candidates[rng.random_range(0..candidates.len())].0;
    }

    let theta_eff = theta * peak;
    let mut total = 0.0;
    for (_, gain) in candidates.iter_mut() {
        total += ((*gain - peak) / theta_eff).exp();
        *gain = total;
    }

    let drawn = rng.random::<f64>() * total;
    candidates
        .iter()
        .find(|(_, cumulative)| drawn < *cumulative)
        .map_or(best, |(community, _)| *community)
}

/// Ties break on the lower community, as in `local_move`. Without it the winner follows
/// whatever order the incident weights happen to be visited in.
fn refine(level: &Level, of: &[usize], gamma: f64, theta: Option<f64>, seed: u64) -> Vec<usize> {
    let mut refined = (0..level.len()).collect::<Vec<_>>();
    let mut sizes = level.size.clone();
    let outer = community_sizes(level, of);
    let mut incident = Incident::new(level.len());
    let mut rng = theta.map(|_| StdRng::seed_from_u64(seed.wrapping_add(SAMPLE_SEED_OFFSET)));
    let mut admitted = Vec::new();

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
        admitted.clear();
        admitted.push((refined[node], 0.0));
        for (candidate, weight) in incident.iter() {
            if candidate == refined[node] || of[candidate] != community {
                continue;
            }
            let gain = weight - gamma * level.size[node] * sizes[candidate];
            if gain > best_gain || (gain == best_gain && candidate < best) {
                best = candidate;
                best_gain = gain;
            }
            if gain >= 0.0 {
                admitted.push((candidate, gain));
            }
        }

        let target = match (theta, rng.as_mut()) {
            (Some(theta), Some(rng)) => draw(&mut admitted, best, best_gain, theta, rng),
            _ => best,
        };
        if target != refined[node] {
            sizes[refined[node]] -= level.size[node];
            sizes[target] += level.size[node];
            refined[node] = target;
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

pub fn leiden(graph: &Graph, gamma: f64, theta: Option<f64>, seed: u64) -> Vec<i32> {
    let mut level = Level::from_graph(graph);
    let mut membership = (0..graph.rows()).collect::<Vec<_>>();
    let mut start: Option<Vec<usize>> = None;
    let mut labels = vec![0i32; graph.rows()];

    for _ in 0..MAX_LEVELS {
        let of = local_move(&level, gamma, seed, start.as_deref());
        for (node, slot) in labels.iter_mut().enumerate() {
            *slot = of[membership[node]] as i32;
        }

        let refined = refine(&level, &of, gamma, theta, seed);
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

pub fn resolutions(graph: &Graph, steps: usize) -> Vec<f64> {
    let weights = Weights::of(graph);
    let nodes = graph.rows().max(1) as f64;
    let mean_degree = 2.0 * weights.total / nodes;
    let largest = (nodes / 2.0).max(2.0);
    let smallest = (nodes / 512.0).max(2.0);
    if steps < 2 || largest <= smallest {
        return vec![mean_degree / largest];
    }
    let ratio = (smallest / largest).powf(1.0 / (steps - 1) as f64);
    (0..steps)
        .map(|step| mean_degree / (largest * ratio.powi(step as i32)))
        .collect()
}
