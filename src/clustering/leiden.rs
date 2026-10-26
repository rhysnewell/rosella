use crate::clustering::graph_partition::{Incident, compact, edge_weight_total, visit_order};
use crate::clustering::layout::{Base, Level};
use crate::embedding::Graph;

const MAX_LEVELS: usize = 20;

fn rank(base: Option<&Base>, node: usize) -> usize {
    base.map_or(node, |base| base.rank(node))
}

fn visits(level: &Level, seed: u64, base: Option<&Base>) -> Vec<usize> {
    base.map_or_else(|| visit_order(level.len(), seed), |base| base.visits(seed))
}

fn community_sizes(level: &Level, of: &[u32]) -> Vec<f64> {
    let mut sizes = vec![0.0; level.len()];
    for (node, community) in of.iter().enumerate() {
        sizes[*community as usize] += level.size[node];
    }
    sizes
}

// Traag, Waltman & van Eck (2019) offer every node an empty community at gain zero, so a node
// whose own community now costs more than it holds can leave it. Communities are 32 bit so the
// lookups a gather scatters over the graph cover half the memory.
fn local_move(
    level: &Level,
    gamma: f64,
    seed: u64,
    start: Option<&[u32]>,
    base: Option<&Base>,
) -> Vec<u32> {
    let everyone = || (0..level.len() as u32).collect::<Vec<_>>();
    let mut of = start.map_or_else(everyone, <[u32]>::to_vec);
    let mut sizes = community_sizes(level, &of);
    let mut members = vec![0usize; level.len()];
    for community in &of {
        members[*community as usize] += 1;
    }
    let mut empty = (0..level.len())
        .rev()
        .filter(|community| members[*community] == 0)
        .collect::<Vec<_>>();
    let order = visits(level, seed, base);
    let mut queued = vec![true; level.len()];
    let mut queue = order
        .iter()
        .copied()
        .collect::<std::collections::VecDeque<_>>();
    let mut incident = Incident::new(level.len());

    while let Some(node) = queue.pop_front() {
        queued[node] = false;
        let current = of[node] as usize;
        level.gather(&mut incident, node, &of);
        sizes[current] -= level.size[node];
        members[current] -= 1;
        if members[current] == 0 {
            empty.push(current);
        }

        let mut best = current;
        let mut best_gain = incident.get(current) - gamma * level.size[node] * sizes[current];
        for (community, weight) in incident.iter() {
            let gain = weight - gamma * level.size[node] * sizes[community];
            if gain > best_gain || (gain == best_gain && rank(base, community) < rank(base, best)) {
                best = community;
                best_gain = gain;
            }
        }
        if let Some(vacant) = empty.last().copied()
            && (0.0 > best_gain || (0.0 == best_gain && rank(base, vacant) < rank(base, best)))
        {
            best = vacant;
        }

        if members[best] == 0 {
            empty.pop();
        }
        sizes[best] += level.size[node];
        members[best] += 1;
        if best != current {
            of[node] = best as u32;
            for (target, _) in level.neighbours(node) {
                if of[target] as usize != best && !queued[target] {
                    queued[target] = true;
                    queue.push_back(target);
                }
            }
        }
    }
    of
}

// Ties break on the lower community, as in `local_move`. Without it the winner follows
// whatever order the incident weights happen to be visited in.
fn refine(level: &Level, of: &[u32], gamma: f64, seed: u64, base: Option<&Base>) -> Vec<u32> {
    let mut refined = (0..level.len() as u32).collect::<Vec<_>>();
    let mut sizes = level.size.clone();
    let outer = community_sizes(level, of);
    let mut incident = Incident::new(level.len());

    for node in visits(level, seed.wrapping_add(1), base) {
        if sizes[refined[node] as usize] != level.size[node] {
            continue;
        }
        let community = of[node];
        let outward = level
            .neighbours(node)
            .filter(|(target, _)| of[*target] == community)
            .map(|(_, weight)| weight)
            .sum::<f64>();
        if outward < gamma * level.size[node] * (outer[community as usize] - level.size[node]) {
            continue;
        }

        let mut best = refined[node] as usize;
        let mut best_gain = 0.0;
        level.gather(&mut incident, node, &refined);
        for (candidate, weight) in incident.iter() {
            if candidate == refined[node] as usize || of[candidate] != community {
                continue;
            }
            let gain = weight - gamma * level.size[node] * sizes[candidate];
            if gain > best_gain || (gain == best_gain && rank(base, candidate) < rank(base, best)) {
                best = candidate;
                best_gain = gain;
            }
        }

        if best != refined[node] as usize {
            sizes[refined[node] as usize] -= level.size[node];
            sizes[best] += level.size[node];
            refined[node] = best as u32;
        }
    }
    refined
}

fn aggregate(level: &Level, refined: &[u32], base: Option<&Base>) -> (Level, Vec<usize>) {
    let nodes = base.map_or_else(|| (0..level.len()).collect(), Base::in_given_order);
    let mut first = vec![usize::MAX; level.len()];
    let mut ids = vec![0usize; level.len()];
    let mut count = 0;
    for node in &nodes {
        let id = &mut first[refined[*node] as usize];
        if *id == usize::MAX {
            *id = count;
            count += 1;
        }
        ids[*node] = *id;
    }

    let mut size = vec![0.0; count];
    let mut members = vec![Vec::new(); count];
    for node in nodes {
        size[ids[node]] += level.size[node];
        members[ids[node]].push(node);
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

// Every rung starts from the same first level, so the rungs share it rather than each holding
// a copy of the whole graph.
pub fn leiden_from(base: &Base, gamma: f64, seed: u64) -> Vec<u32> {
    let mut aggregated: Option<Level> = None;
    let mut membership = base.in_given_order();
    let mut start: Option<Vec<u32>> = None;
    let mut labels = vec![0u32; membership.len()];

    for _ in 0..MAX_LEVELS {
        let (level, laid) = match &aggregated {
            Some(level) => (level, None),
            None => (&base.level, Some(base)),
        };
        let of = local_move(level, gamma, seed, start.as_deref(), laid);
        for (node, slot) in labels.iter_mut().enumerate() {
            *slot = of[membership[node]];
        }

        let refined = refine(level, &of, gamma, seed, laid);
        let (next, ids) = aggregate(level, &refined, laid);
        if next.len() == level.len() {
            break;
        }

        let mut inherited = vec![0u32; next.len()];
        for (node, id) in ids.iter().enumerate() {
            inherited[*id] = of[node];
        }
        start = Some(compact(&inherited));

        for slot in &mut membership {
            *slot = ids[*slot];
        }
        aggregated = Some(next);
    }

    compact(&labels)
}

pub fn resolutions(graph: &Graph, sizes: &[f64], steps: usize) -> Vec<f64> {
    let edges = edge_weight_total(graph);
    let total = sizes.iter().sum::<f64>().max(1.0);
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
