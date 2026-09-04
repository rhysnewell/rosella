use std::collections::HashMap;

use statrs::function::gamma::ln_gamma;

fn xlogx(value: f64) -> f64 {
    if value > 0.0 { value * value.ln() } else { 0.0 }
}

fn pair_term(a: usize, b: usize, weight: f64) -> f64 {
    if a == b {
        -0.5 * xlogx(weight)
    } else {
        -xlogx(weight)
    }
}

/// The pairs of the block matrix a change to `first` and `second` can reach, each once.
fn pair_list(touched: &[usize], first: usize, second: usize) -> Vec<(usize, usize)> {
    let mut pairs = Vec::with_capacity(touched.len() * 2);
    for target in touched {
        for source in [first, second] {
            if *target < source && (*target == first || *target == second) {
                continue;
            }
            pairs.push((source, *target));
        }
    }
    pairs
}

/// A node's weight into each block it reaches, gathered once per visit and sorted by block
/// so `move_deltas` can look a block up without scanning.
pub type Incidence = [(usize, f64)];

/// Weighted edge counts between blocks, the block degrees and the block sizes. Every term of
/// the description length reads only these, so a change costs the blocks it touches and no more.
pub struct BlockState {
    pub of: Vec<usize>,
    rows: Vec<HashMap<usize, f64>>,
    degree: Vec<f64>,
    count: Vec<usize>,
    blocks: usize,
    nodes: usize,
    edges: f64,
}

impl BlockState {
    pub fn new(
        neighbours: &[Vec<(usize, f64)>],
        node_degree: &[f64],
        of: Vec<usize>,
        block_count: usize,
    ) -> Self {
        let nodes = neighbours.len();
        let mut state = Self {
            of,
            rows: vec![HashMap::new(); block_count],
            degree: vec![0.0; block_count],
            count: vec![0; block_count],
            blocks: 0,
            nodes,
            edges: node_degree.iter().sum::<f64>() / 2.0,
        };
        for node in 0..nodes {
            let block = state.of[node];
            state.degree[block] += node_degree[node];
            state.count[block] += 1;
            for (target, weight) in &neighbours[node] {
                *state.rows[block].entry(state.of[*target]).or_insert(0.0) += weight;
            }
        }
        state.blocks = state.count.iter().filter(|count| **count > 0).count();
        state
    }

    pub fn blocks(&self) -> usize {
        self.blocks
    }

    pub fn count(&self, block: usize) -> usize {
        self.count[block]
    }

    pub fn neighbouring_blocks(&self, block: usize) -> impl Iterator<Item = usize> + '_ {
        self.rows[block].keys().copied()
    }

    fn pair(&self, a: usize, b: usize) -> f64 {
        self.rows[a].get(&b).copied().unwrap_or(0.0)
    }

    fn add(&mut self, a: usize, b: usize, delta: f64) {
        if delta == 0.0 {
            return;
        }
        for (row, column) in if a == b {
            vec![(a, b)]
        } else {
            vec![(a, b), (b, a)]
        } {
            let value = self.rows[row].entry(column).or_insert(0.0);
            *value += delta;
            if *value <= 0.0 {
                self.rows[row].remove(&column);
            }
        }
    }

    /// Peixoto (2014). The block matrix and degree sequence priors are what let the block count
    /// come out of the fit rather than out of a resolution parameter.
    fn global(&self, blocks: usize) -> f64 {
        if blocks == 0 {
            return 0.0;
        }
        let nodes = self.nodes as f64;
        let pairs = blocks as f64 * (blocks as f64 + 1.0) / 2.0;
        let matrix = ln_gamma(self.edges + pairs) - ln_gamma(self.edges + 1.0) - ln_gamma(pairs);
        let counts =
            ln_gamma(nodes) - ln_gamma(blocks as f64) - ln_gamma(nodes - blocks as f64 + 1.0);
        matrix + counts + nodes.ln()
    }

    fn block_term(&self, degree: f64, count: usize) -> f64 {
        if count == 0 {
            return 0.0;
        }
        let count = count as f64;
        xlogx(degree) - ln_gamma(count + 1.0) + ln_gamma(degree + count)
            - ln_gamma(degree + 1.0)
            - ln_gamma(count)
    }

    pub fn description_length(&self) -> f64 {
        let mut total = self.global(self.blocks);
        for block in 0..self.rows.len() {
            total += self.block_term(self.degree[block], self.count[block]);
            for (other, weight) in &self.rows[block] {
                if *other >= block {
                    total += pair_term(block, *other, *weight);
                }
            }
        }
        total
    }

    fn touched(
        &self,
        first: usize,
        second: usize,
        reach: impl Iterator<Item = usize>,
    ) -> Vec<usize> {
        let mut touched = vec![first, second];
        touched.extend(reach);
        touched.sort_unstable();
        touched.dedup();
        touched
    }

    fn cost(&self, pairs: &[(usize, usize)], deltas: &[f64], ends: [(usize, f64, i64); 2]) -> f64 {
        let mut before = 0.0;
        let mut after = 0.0;
        for ((source, target), delta) in pairs.iter().zip(deltas) {
            let current = self.pair(*source, *target);
            before += pair_term(*source, *target, current);
            after += pair_term(*source, *target, current + delta);
        }

        let mut blocks = self.blocks;
        for (block, degree_delta, count_delta) in ends {
            before += self.block_term(self.degree[block], self.count[block]);
            let count = (self.count[block] as i64 + count_delta) as usize;
            after += self.block_term(self.degree[block] + degree_delta, count);
            if self.count[block] > 0 && count == 0 {
                blocks -= 1;
            }
            if self.count[block] == 0 && count > 0 {
                blocks += 1;
            }
        }

        after + self.global(blocks) - before - self.global(self.blocks)
    }

    fn move_deltas(
        &self,
        pairs: &[(usize, usize)],
        from: usize,
        to: usize,
        incidence: &Incidence,
    ) -> Vec<f64> {
        let weight_of = |block: usize| {
            incidence
                .binary_search_by_key(&block, |(other, _)| *other)
                .map_or(0.0, |at| incidence[at].1)
        };
        let (into_from, into_to) = (weight_of(from), weight_of(to));
        pairs
            .iter()
            .map(|(source, target)| match (*source, *target) {
                (a, b) if a == from && b == from => -2.0 * into_from,
                (a, b) if a == to && b == to => 2.0 * into_to,
                (a, b) if (a == from && b == to) || (a == to && b == from) => into_from - into_to,
                (a, b) if a == from => -weight_of(b),
                (a, b) if a == to => weight_of(b),
                _ => 0.0,
            })
            .collect()
    }

    fn merge_deltas(&self, pairs: &[(usize, usize)], keep: usize, gone: usize) -> Vec<f64> {
        let between = self.pair(keep, gone);
        let inside = self.pair(gone, gone);
        pairs
            .iter()
            .map(|(source, target)| match (*source, *target) {
                (a, b) if a == keep && b == keep => 2.0 * between + inside,
                (a, b) if a == gone && b == gone => -inside,
                (a, b) if (a == keep && b == gone) || (a == gone && b == keep) => -between,
                (a, b) if a == keep => self.pair(gone, b),
                (a, b) if a == gone => -self.pair(gone, b),
                _ => 0.0,
            })
            .collect()
    }

    pub fn move_cost(&self, node: usize, to: usize, incidence: &Incidence, degree: f64) -> f64 {
        let from = self.of[node];
        if from == to {
            return 0.0;
        }
        let touched = self.touched(from, to, incidence.iter().map(|(block, _)| *block));
        let pairs = pair_list(&touched, from, to);
        let deltas = self.move_deltas(&pairs, from, to, incidence);
        self.cost(&pairs, &deltas, [(from, -degree, -1), (to, degree, 1)])
    }

    pub fn apply_move(&mut self, node: usize, to: usize, incidence: &Incidence, degree: f64) {
        let from = self.of[node];
        if from == to {
            return;
        }
        let touched = self.touched(from, to, incidence.iter().map(|(block, _)| *block));
        let pairs = pair_list(&touched, from, to);
        let deltas = self.move_deltas(&pairs, from, to, incidence);
        for ((source, target), delta) in pairs.iter().zip(&deltas) {
            self.add(*source, *target, *delta);
        }
        if self.count[from] == 1 {
            self.blocks -= 1;
        }
        if self.count[to] == 0 {
            self.blocks += 1;
        }
        self.degree[from] -= degree;
        self.degree[to] += degree;
        self.count[from] -= 1;
        self.count[to] += 1;
        self.of[node] = to;
    }

    fn merge_pairs(&self, keep: usize, gone: usize) -> Vec<(usize, usize)> {
        let reach = self.rows[keep]
            .keys()
            .chain(self.rows[gone].keys())
            .copied()
            .collect::<Vec<_>>();
        let touched = self.touched(keep, gone, reach.into_iter());
        pair_list(&touched, keep, gone)
    }

    pub fn merge_cost(&self, keep: usize, gone: usize) -> f64 {
        if keep == gone || self.count[gone] == 0 {
            return f64::INFINITY;
        }
        let pairs = self.merge_pairs(keep, gone);
        let deltas = self.merge_deltas(&pairs, keep, gone);
        let moved = self.count[gone] as i64;
        self.cost(
            &pairs,
            &deltas,
            [
                (gone, -self.degree[gone], -moved),
                (keep, self.degree[gone], moved),
            ],
        )
    }

    pub fn apply_merge(&mut self, keep: usize, gone: usize) {
        if keep == gone || self.count[gone] == 0 {
            return;
        }
        let pairs = self.merge_pairs(keep, gone);
        let deltas = self.merge_deltas(&pairs, keep, gone);
        for ((source, target), delta) in pairs.iter().zip(&deltas) {
            self.add(*source, *target, *delta);
        }
        self.degree[keep] += self.degree[gone];
        self.count[keep] += self.count[gone];
        self.degree[gone] = 0.0;
        self.count[gone] = 0;
        self.blocks -= 1;
    }

    pub fn capacity(&self) -> usize {
        self.rows.len()
    }

    pub fn live_blocks(&self) -> Vec<usize> {
        (0..self.capacity())
            .filter(|block| self.count[*block] > 0)
            .collect()
    }

    /// Merges leave `of` stale on purpose: rewriting it per merge is a pass over every node,
    /// and a round does up to half the blocks.
    pub fn relabel(&mut self, redirect: &[usize]) {
        for block in &mut self.of {
            let mut target = *block;
            while redirect[target] != target {
                target = redirect[target];
            }
            *block = target;
        }
    }
}
