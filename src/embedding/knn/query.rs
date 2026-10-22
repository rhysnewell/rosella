use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use super::{KnnGraph, MAX_ITERATIONS, Metric, NeighbourList, graph_of};
use crate::embedding::metrics::prepared::{LANES, PreparedAggregate};

// Each query's k nearest points of a base whose own graph is already built, found by walking
// that graph from random starts. Queries never meet, so each is searched alone and in parallel.
pub fn nearest_in<M: Metric>(
    base: &KnnGraph,
    queries: usize,
    k: usize,
    max_candidates: usize,
    seed: u64,
    metric: M,
) -> KnnGraph {
    let n = base.n_points();
    let k = k.min(n).max(1);
    if n == 0 {
        return KnnGraph {
            indices: Array2::from_elem((queries, k), u32::MAX),
            dists: Array2::from_elem((queries, k), f32::INFINITY),
        };
    }
    let lists = (0..queries)
        .into_par_iter()
        .map(|query| {
            let mut rng = StdRng::seed_from_u64(
                seed ^ (query as u64).wrapping_mul(crate::defaults::SEED_STRIDE),
            );
            let mut list = NeighbourList::new(k);
            let mut seen = Seen::default();
            let (mut unseen, mut distances) = (Vec::new(), Vec::new());
            for _ in 0..k {
                let start = rng.random_range(0..n);
                if seen.insert(start as u32) {
                    list.push(metric.distance(query, start), start as u32);
                }
            }
            for _ in 0..MAX_ITERATIONS {
                let fresh = (0..k)
                    .filter(|slot| list.is_new[*slot])
                    .take(max_candidates)
                    .collect::<Vec<_>>();
                if fresh.is_empty() {
                    break;
                }
                let from = fresh
                    .into_iter()
                    .map(|slot| {
                        list.is_new[slot] = false;
                        list.indices[slot]
                    })
                    .collect::<Vec<_>>();
                for point in from {
                    unseen.clear();
                    distances.clear();
                    unseen.extend(
                        base.indices
                            .row(point as usize)
                            .iter()
                            .take(max_candidates)
                            .filter(|next| **next != u32::MAX && seen.insert(**next)),
                    );
                    let bound = list.worst();
                    metric.within_many(query, &unseen, |_| bound, &mut distances);
                    for (next, distance) in unseen.iter().zip(&distances) {
                        list.push(*distance, *next);
                    }
                }
            }
            list
        })
        .collect::<Vec<_>>();
    graph_of(&lists, k)
}

// Exact where the descent is approximate, for the few rows a decision rests on.
pub fn nearest_exact(
    base: usize,
    queries: &[usize],
    k: usize,
    metric: &PreparedAggregate,
) -> KnnGraph {
    let k = k.min(base.saturating_sub(1)).max(1);
    let lists = queries
        .par_chunks(LANES)
        .flat_map_iter(|queries| {
            let block = metric.block(queries);
            let mut lists = queries
                .iter()
                .map(|_| NeighbourList::new(k))
                .collect::<Vec<_>>();
            let (mut bounds, mut distances) = ([0.0; LANES], [0.0; LANES]);
            for other in 0..base {
                for (bound, list) in bounds.iter_mut().zip(&lists) {
                    *bound = list.worst();
                }
                metric.within_block(&block, other, &bounds, &mut distances);
                for ((query, list), distance) in queries.iter().zip(&mut lists).zip(distances) {
                    if *query != other {
                        list.push(distance, other as u32);
                    }
                }
            }
            lists
        })
        .collect::<Vec<_>>();
    graph_of(&lists, k)
}

// A walk inserts every point it meets, so the set hashes on the search's hottest path.
type Seen = std::collections::HashSet<u32, std::hash::BuildHasherDefault<Scatter>>;

#[derive(Default)]
struct Scatter(u64);

impl std::hash::Hasher for Scatter {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for byte in bytes {
            self.0 = (self.0.rotate_left(8) ^ u64::from(*byte)).wrapping_mul(SCATTER);
        }
    }

    fn write_u32(&mut self, value: u32) {
        self.0 = u64::from(value).wrapping_mul(SCATTER);
    }
}

const SCATTER: u64 = 0x9E37_79B9_7F4A_7C15;
