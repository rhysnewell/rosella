use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;
use std::sync::Mutex;

const MAX_CANDIDATES: usize = 32;
const MAX_ITERATIONS: usize = 20;
const CONVERGENCE_FRACTION: f64 = 0.001;
const ROW_SEED_STRIDE: u64 = 0x9E37_79B9_7F4A_7C15;

pub struct KnnGraph {
    pub indices: Array2<u32>,
    pub dists: Array2<f32>,
}

impl KnnGraph {
    pub fn n_points(&self) -> usize {
        self.indices.nrows()
    }

    /// Points whose every neighbour slot stayed empty.
    pub fn disconnected(&self) -> Vec<usize> {
        (0..self.n_points())
            .filter(|row| {
                self.indices
                    .row(*row)
                    .iter()
                    .all(|index| *index == u32::MAX)
            })
            .collect()
    }
}

/// A bounded set of the k closest neighbours seen so far, sorted by distance and then
/// index. Insertion is order independent, which is what makes the whole build
/// reproducible under `rayon`.
struct NeighbourList {
    dists: Vec<f64>,
    indices: Vec<u32>,
    is_new: Vec<bool>,
}

impl NeighbourList {
    fn new(k: usize) -> Self {
        Self {
            dists: vec![f64::INFINITY; k],
            indices: vec![u32::MAX; k],
            is_new: vec![false; k],
        }
    }

    /// Total order on (distance, index). The index tiebreak is what keeps insertion
    /// order independent.
    fn sorts_before(&self, position: usize, distance: f64, index: u32) -> bool {
        self.dists[position] < distance
            || (self.dists[position] == distance && self.indices[position] < index)
    }

    fn push(&mut self, distance: f64, index: u32) -> bool {
        let k = self.dists.len();
        if self.sorts_before(k - 1, distance, index) {
            return false;
        }
        if self.indices.contains(&index) {
            return false;
        }

        let mut position = k - 1;
        while position > 0 && !self.sorts_before(position - 1, distance, index) {
            self.dists[position] = self.dists[position - 1];
            self.indices[position] = self.indices[position - 1];
            self.is_new[position] = self.is_new[position - 1];
            position -= 1;
        }
        self.dists[position] = distance;
        self.indices[position] = index;
        self.is_new[position] = true;
        true
    }
}

/// Nearest neighbour descent. Deterministic for a given seed and k, whatever the thread
/// count, because the only randomness is the seeded initial sample and every later step
/// is order independent.
///
/// The metric takes row indices rather than rows, because anything read off the contig
/// itself, its length now and its markers later, needs to know which contig it is looking at.
pub fn build_knn<M>(n: usize, k: usize, seed: u64, metric: M) -> KnnGraph
where
    M: Fn(usize, usize) -> f64 + Sync,
{
    let k = k.min(n.saturating_sub(1)).max(1);

    let neighbours = (0..n)
        .map(|_| Mutex::new(NeighbourList::new(k)))
        .collect::<Vec<_>>();

    neighbours.par_iter().enumerate().for_each(|(i, list)| {
        let mut rng = StdRng::seed_from_u64(seed ^ (i as u64).wrapping_mul(ROW_SEED_STRIDE));
        let mut list = list.lock().unwrap();
        for _ in 0..k {
            let j = rng.random_range(0..n);
            if j != i {
                list.push(metric(i, j), j as u32);
            }
        }
    });

    for _ in 0..MAX_ITERATIONS {
        let (new_candidates, old_candidates) = build_candidates(&neighbours, n);

        let updates: usize = (0..n)
            .into_par_iter()
            .map(|i| join(&metric, &neighbours, &new_candidates[i], &old_candidates[i]))
            .sum();

        if updates as f64 <= CONVERGENCE_FRACTION * k as f64 * n as f64 {
            break;
        }
    }

    let mut indices = Array2::from_elem((n, k), u32::MAX);
    let mut dists = Array2::from_elem((n, k), f32::INFINITY);
    for (i, list) in neighbours.iter().enumerate() {
        let list = list.lock().unwrap();
        for j in 0..k {
            indices[[i, j]] = list.indices[j];
            dists[[i, j]] = list.dists[j] as f32;
        }
    }

    KnnGraph { indices, dists }
}

/// Forward and reverse neighbour lists, split by whether the edge is new since the last
/// pass. Built in index order so the caps fall the same way every run.
fn build_candidates(
    neighbours: &[Mutex<NeighbourList>],
    n: usize,
) -> (Vec<Vec<u32>>, Vec<Vec<u32>>) {
    let mut new_candidates = vec![Vec::new(); n];
    let mut old_candidates = vec![Vec::new(); n];

    for i in 0..n {
        let mut list = neighbours[i].lock().unwrap();
        for slot in 0..list.indices.len() {
            let neighbour = list.indices[slot];
            if neighbour == u32::MAX {
                continue;
            }
            let target = if list.is_new[slot] {
                &mut new_candidates
            } else {
                &mut old_candidates
            };
            if target[i].len() < MAX_CANDIDATES {
                target[i].push(neighbour);
                list.is_new[slot] = false;
            }
            if target[neighbour as usize].len() < MAX_CANDIDATES {
                target[neighbour as usize].push(i as u32);
            }
        }
    }

    (new_candidates, old_candidates)
}

fn join<M>(
    metric: &M,
    neighbours: &[Mutex<NeighbourList>],
    new_candidates: &[u32],
    old_candidates: &[u32],
) -> usize
where
    M: Fn(usize, usize) -> f64 + Sync,
{
    let mut updates = 0;
    for (position, a) in new_candidates.iter().enumerate() {
        let pairs = new_candidates[position + 1..]
            .iter()
            .chain(old_candidates.iter());
        for b in pairs {
            if a == b {
                continue;
            }
            let distance = metric(*a as usize, *b as usize);
            if neighbours[*a as usize].lock().unwrap().push(distance, *b) {
                updates += 1;
            }
            if neighbours[*b as usize].lock().unwrap().push(distance, *a) {
                updates += 1;
            }
        }
    }
    updates
}

/// Exact k nearest neighbours. Quadratic, so it exists to check `build_knn` rather than
/// to run on real assemblies.
pub fn brute_force_knn<M>(n: usize, k: usize, metric: M) -> KnnGraph
where
    M: Fn(usize, usize) -> f64 + Sync,
{
    let k = k.min(n.saturating_sub(1)).max(1);

    let mut indices = Array2::from_elem((n, k), u32::MAX);
    let mut dists = Array2::from_elem((n, k), f32::INFINITY);

    let rows_of_neighbours = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut list = NeighbourList::new(k);
            for j in 0..n {
                if i != j {
                    list.push(metric(i, j), j as u32);
                }
            }
            list
        })
        .collect::<Vec<_>>();

    for (i, list) in rows_of_neighbours.iter().enumerate() {
        for j in 0..k {
            indices[[i, j]] = list.indices[j];
            dists[[i, j]] = list.dists[j] as f32;
        }
    }

    KnnGraph { indices, dists }
}
