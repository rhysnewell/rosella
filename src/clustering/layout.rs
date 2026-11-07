use crate::clustering::graph_partition::{Incident, visit_order};
use crate::embedding::{Graph, row_of};

// One flat run of edges, so a gather walks memory in order instead of chasing a row per node.
pub(super) struct Level {
    offsets: Vec<usize>,
    targets: Vec<u32>,
    weights: Vec<f64>,
    pub(super) size: Vec<f64>,
}

impl Level {
    pub(super) fn empty(size: Vec<f64>) -> Self {
        Self {
            offsets: vec![0],
            targets: Vec::new(),
            weights: Vec::new(),
            size,
        }
    }

    pub(super) fn close(&mut self, edges: impl Iterator<Item = (usize, f64)>) {
        for (target, weight) in edges {
            self.targets.push(target as u32);
            self.weights.push(weight);
        }
        self.offsets.push(self.targets.len());
    }

    pub(super) fn neighbours(&self, node: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let span = self.offsets[node]..self.offsets[node + 1];
        self.targets[span.clone()]
            .iter()
            .zip(&self.weights[span])
            .map(|(target, weight)| (*target as usize, *weight))
    }

    pub(super) fn len(&self) -> usize {
        self.size.len()
    }

    pub(super) fn gather(&self, incident: &mut Incident, node: usize, of: &[u32]) {
        incident.gather(
            self.neighbours(node)
                .map(|(target, weight)| (of[target] as usize, weight)),
        );
    }
}

// A base level numbered for locality still visits and breaks ties by the given numbering, so it
// labels every node as the given numbering would.
pub struct Base {
    pub(super) level: Level,
    original: Vec<u32>,
    placed: Vec<u32>,
    looped: Vec<bool>,
}

impl Base {
    pub(super) fn rank(&self, node: usize) -> usize {
        self.original[node] as usize
    }

    pub(super) fn visits(&self, seed: u64) -> Vec<usize> {
        visit_order(self.placed.len(), seed)
            .into_iter()
            .map(|node| self.placed[node] as usize)
            .collect()
    }

    pub(super) fn in_given_order(&self) -> Vec<usize> {
        self.placed.iter().map(|at| *at as usize).collect()
    }

    // A node on a self loop reads itself, so its own move marks it stale as the graph's rows do.
    pub(super) fn readers(&self) -> (Vec<usize>, Vec<u32>) {
        let nodes = self.level.len();
        let mut starts = vec![0usize; nodes + 1];
        for node in 0..nodes {
            for (target, _) in self.level.neighbours(node) {
                starts[target + 1] += 1;
            }
            starts[node + 1] += usize::from(self.looped[node]);
        }
        for at in 0..nodes {
            starts[at + 1] += starts[at];
        }
        let mut next = starts.clone();
        let mut readers = vec![0u32; starts[nodes]];
        for node in 0..nodes {
            let targets = self.level.neighbours(node).map(|(target, _)| target);
            for target in targets.chain(self.looped[node].then_some(node)) {
                readers[next[target]] = node as u32;
                next[target] += 1;
            }
        }
        (starts, readers)
    }
}

// The walk gives a node's neighbours ids near its own, so the community lookups a gather scatters
// land in a few cache lines rather than across the graph.
pub fn base_level(graph: &Graph, sizes: &[f64]) -> Base {
    base_level_in(graph, sizes, breadth_first(graph))
}

pub fn base_level_in(graph: &Graph, sizes: &[f64], original: Vec<u32>) -> Base {
    let mut placed = vec![0u32; original.len()];
    for (at, node) in original.iter().enumerate() {
        placed[*node as usize] = at as u32;
    }
    let mut level = Level::empty(original.iter().map(|node| sizes[*node as usize]).collect());
    let mut looped = vec![false; original.len()];
    for (at, node) in original.iter().enumerate() {
        let row = *node as usize;
        let (targets, weights) = row_of(graph, row);
        looped[at] = targets.contains(node);
        level.close(
            targets
                .iter()
                .zip(weights)
                .filter(|(target, _)| **target as usize != row)
                .map(|(target, weight)| (placed[*target as usize] as usize, *weight as f64)),
        );
    }
    Base {
        level,
        original,
        placed,
        looped,
    }
}

fn breadth_first(graph: &Graph) -> Vec<u32> {
    let mut seen = vec![false; graph.rows()];
    let mut order = Vec::with_capacity(graph.rows());
    for root in 0..graph.rows() {
        if std::mem::replace(&mut seen[root], true) {
            continue;
        }
        let mut next = order.len();
        order.push(root as u32);
        while let Some(node) = order.get(next).copied() {
            next += 1;
            for target in row_of(graph, node as usize).0 {
                if !std::mem::replace(&mut seen[*target as usize], true) {
                    order.push(*target);
                }
            }
        }
    }
    order
}
