use std::sync::{
    Mutex,
    atomic::{AtomicU64, Ordering},
};

use log::debug;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use super::{KnnGraph, MAX_ITERATIONS, Metric, NeighbourList, candidates::Candidates, graph_of};

// A list's worst only falls, so an offer past a stale read of it is past the list too and is
// refused without the lock. What the lists keep is unchanged.
struct Lists {
    rows: Vec<Mutex<NeighbourList>>,
    worst: Vec<AtomicU64>,
}

impl Lists {
    fn new(n: usize, k: usize) -> Self {
        Self {
            rows: (0..n).map(|_| Mutex::new(NeighbourList::new(k))).collect(),
            worst: (0..n)
                .map(|_| AtomicU64::new(f64::INFINITY.to_bits()))
                .collect(),
        }
    }

    fn bound(&self, row: usize) -> f64 {
        f64::from_bits(self.worst[row].load(Ordering::Relaxed))
    }

    fn offer(&self, row: usize, distance: f64, index: u32) {
        if distance > self.bound(row) {
            return;
        }
        let mut list = self.rows[row].lock().unwrap();
        list.push(distance, index);
        self.worst[row].store(list.worst().to_bits(), Ordering::Relaxed);
    }

    // A join offers its own row a whole bucket at once, so that row is locked once, not per pair.
    fn offer_many(&self, row: usize, offers: impl Iterator<Item = (f64, u32)>) {
        let bound = self.bound(row);
        let mut offers = offers.filter(|(distance, _)| *distance <= bound).peekable();
        if offers.peek().is_none() {
            return;
        }
        let mut list = self.rows[row].lock().unwrap();
        for (distance, index) in offers {
            list.push(distance, index);
        }
        self.worst[row].store(list.worst().to_bits(), Ordering::Relaxed);
    }
}

// Deterministic for a given seed and k, whatever the thread count, because the lists keep
// the k smallest under a total order and the stop rule reads them rather than the pushes.
pub fn build_knn_with<M: Metric>(
    n: usize,
    k: usize,
    max_candidates: usize,
    seed: u64,
    metric: M,
) -> KnnGraph {
    let k = k.min(n.saturating_sub(1)).max(1);
    let lists = Lists::new(n, k);
    (0..n).into_par_iter().for_each(|i| {
        let mut list = lists.rows[i].lock().unwrap();
        fill_at_random(&mut list, i, n, k, seed, &metric);
        lists.worst[i].store(list.worst().to_bits(), Ordering::Relaxed);
    });
    descend(lists, n, k, max_candidates, seed, metric)
}

fn fill_at_random<M: Metric>(
    list: &mut NeighbourList,
    i: usize,
    n: usize,
    k: usize,
    seed: u64,
    metric: &M,
) {
    let mut rng =
        StdRng::seed_from_u64(seed ^ (i as u64).wrapping_mul(crate::defaults::SEED_STRIDE));
    for _ in 0..k {
        let j = rng.random_range(0..n);
        if j != i {
            list.push(metric.distance(i, j), j as u32);
        }
    }
}

fn descend<M: Metric>(
    lists: Lists,
    n: usize,
    k: usize,
    max_candidates: usize,
    seed: u64,
    metric: M,
) -> KnnGraph {
    let progress = crate::progress::counted(
        crate::progress::Stage::NearestNeighbours,
        MAX_ITERATIONS as u64,
    );
    let mut candidates = Candidates::new(n, max_candidates.max(1));
    for round in 0..MAX_ITERATIONS {
        let drawn = seed ^ (round as u64 + 1).wrapping_mul(crate::defaults::SEED_STRIDE);
        candidates.build(&lists.rows, drawn);

        (0..n)
            .into_par_iter()
            .for_each(|i| join(&metric, &lists, candidates.new_of(i), candidates.old_of(i)));

        let taken = taken_slots(&lists.rows);
        debug!("Descent round {round}: {taken} slots");
        progress.inc(1);
        progress.set_message(format!("{taken} neighbours taken"));
        let fraction =
            crate::tuning::arm("ROSELLA_CONVERGENCE", crate::tuning::CONVERGENCE_FRACTION);
        if taken as f64 <= fraction * k as f64 * n as f64 {
            break;
        }
    }
    crate::progress::done(&progress);

    let rows = lists
        .rows
        .into_iter()
        .map(|row| row.into_inner().unwrap())
        .collect::<Vec<_>>();
    graph_of(&rows, k)
}

fn join<M: Metric>(metric: &M, lists: &Lists, new_candidates: &[u32], old_candidates: &[u32]) {
    let mut others = Vec::with_capacity(new_candidates.len() + old_candidates.len());
    let mut distances = Vec::with_capacity(others.capacity());
    for (position, a) in new_candidates.iter().enumerate() {
        others.clear();
        distances.clear();
        others.extend(
            new_candidates[position + 1..]
                .iter()
                .chain(old_candidates)
                .filter(|b| *b != a),
        );
        let a = *a as usize;
        metric.within_many(
            a,
            &others,
            |b| lists.bound(a).max(lists.bound(b)),
            &mut distances,
        );
        for (b, distance) in others.iter().zip(&distances) {
            lists.offer(*b as usize, *distance, a as u32);
        }
        lists.offer_many(a, distances.iter().copied().zip(others.iter().copied()));
    }
}

// What the lists kept, not what the pushes won. A push that takes a slot and is displaced
// later in the same pass counts to the pushes, and the pushes are the half that races.
fn taken_slots(neighbours: &[Mutex<NeighbourList>]) -> usize {
    neighbours
        .par_iter()
        .map(|list| {
            list.lock()
                .unwrap()
                .fresh
                .iter()
                .filter(|fresh| **fresh)
                .count()
        })
        .sum()
}
