use std::sync::{
    Mutex,
    atomic::{AtomicU64, Ordering},
};

use log::debug;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use super::{KnnGraph, MAX_ITERATIONS, Metric, NeighbourList, graph_of};

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
        build_candidates(&lists.rows, n, &mut candidates, drawn);

        (0..n)
            .into_par_iter()
            .for_each(|i| join(&metric, &lists, candidates.new_of(i), candidates.old_of(i)));

        let taken = taken_slots(&lists.rows);
        debug!("Descent round {round}: {taken} slots");
        progress.inc(1);
        progress.set_message(format!("{taken} neighbours taken"));
        if taken as f64 <= crate::tuning::CONVERGENCE_FRACTION * k as f64 * n as f64 {
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

// Forward and reverse neighbour lists, split by whether the edge is new since the last
// pass. Every bucket is capped, so one flat allocation of that stride holds the whole pass.
struct Candidates {
    stride: usize,
    new: Vec<u32>,
    new_rank: Vec<u32>,
    new_len: Vec<u32>,
    old: Vec<u32>,
    old_rank: Vec<u32>,
    old_len: Vec<u32>,
}

// A bucket keeps the lowest ranks offered, as a max heap on rank, so what it holds does not
// hang on the order the rows offered in.
fn offer_candidate(
    values: &mut [u32],
    ranks: &mut [u32],
    lengths: &mut [u32],
    stride: usize,
    row: usize,
    value: u32,
    rank: u32,
) {
    let values = &mut values[row * stride..(row + 1) * stride];
    let ranks = &mut ranks[row * stride..(row + 1) * stride];
    let length = lengths[row] as usize;
    if length < stride {
        lengths[row] += 1;
        values[length] = value;
        ranks[length] = rank;
        let mut at = length;
        while at > 0 && ranks[(at - 1) / 2] < ranks[at] {
            values.swap(at, (at - 1) / 2);
            ranks.swap(at, (at - 1) / 2);
            at = (at - 1) / 2;
        }
    } else if rank < ranks[0] {
        values[0] = value;
        ranks[0] = rank;
        let mut at = 0;
        loop {
            let largest = [2 * at + 1, 2 * at + 2]
                .into_iter()
                .filter(|child| *child < stride)
                .fold(at, |best, child| match ranks[child] > ranks[best] {
                    true => child,
                    false => best,
                });
            if largest == at {
                break;
            }
            values.swap(at, largest);
            ranks.swap(at, largest);
            at = largest;
        }
    }
}

impl Candidates {
    fn new(n: usize, stride: usize) -> Self {
        Self {
            stride,
            new: vec![u32::MAX; n * stride],
            new_rank: vec![u32::MAX; n * stride],
            new_len: vec![0; n],
            old: vec![u32::MAX; n * stride],
            old_rank: vec![u32::MAX; n * stride],
            old_len: vec![0; n],
        }
    }

    fn offer(&mut self, fresh: bool, row: usize, value: u32, rank: u32) {
        let (values, ranks, lengths) = match fresh {
            true => (&mut self.new, &mut self.new_rank, &mut self.new_len),
            false => (&mut self.old, &mut self.old_rank, &mut self.old_len),
        };
        offer_candidate(values, ranks, lengths, self.stride, row, value, rank);
    }

    fn clear(&mut self) {
        self.new_len.fill(0);
        self.old_len.fill(0);
    }

    fn new_of(&self, row: usize) -> &[u32] {
        &self.new[row * self.stride..row * self.stride + self.new_len[row] as usize]
    }

    fn old_of(&self, row: usize) -> &[u32] {
        &self.old[row * self.stride..row * self.stride + self.old_len[row] as usize]
    }

    fn dedup(&mut self) {
        self.new
            .par_chunks_mut(self.stride)
            .zip(self.new_len.par_iter_mut())
            .zip(
                self.old
                    .par_chunks_mut(self.stride)
                    .zip(self.old_len.par_iter_mut()),
            )
            .for_each(|((new, new_len), (old, old_len))| {
                *new_len = unique(&mut new[..*new_len as usize], &[]);
                *old_len = unique(&mut old[..*old_len as usize], &new[..*new_len as usize]);
            });
    }
}

// Dong, Charikar & Li (2011) sample each bucket at random. Filled first come by row, late rows
// starved and the search converged as far as the contig order let it.
fn build_candidates(
    neighbours: &[Mutex<NeighbourList>],
    n: usize,
    candidates: &mut Candidates,
    seed: u64,
) {
    candidates.clear();
    for (i, entry) in neighbours.iter().enumerate().take(n) {
        let mut list = entry.lock().unwrap();
        list.fresh.fill(false);
        for slot in 0..list.indices.len() {
            let neighbour = list.indices[slot];
            if neighbour == u32::MAX {
                continue;
            }
            let rank =
                (crate::digest::mix(seed ^ ((u64::from(neighbour) << 32) | i as u64)) >> 32) as u32;
            candidates.offer(list.is_new[slot], i, neighbour, rank);
            candidates.offer(list.is_new[slot], neighbour as usize, i as u32, rank);
        }
    }
    candidates.dedup();
    neighbours.par_iter().enumerate().for_each(|(i, entry)| {
        let mut list = entry.lock().unwrap();
        let drawn = candidates.new_of(i);
        for slot in 0..list.indices.len() {
            if list.is_new[slot] && drawn.binary_search(&list.indices[slot]).is_ok() {
                list.is_new[slot] = false;
            }
        }
    });
}

// A mutual edge reaches a bucket from both ends and a contig can sit in both buckets, so the join
// would measure one pair twice. The caps are already spent, so the lists come out the same.
fn unique(values: &mut [u32], elsewhere: &[u32]) -> u32 {
    values.sort_unstable();
    let mut kept = 0;
    for at in 0..values.len() {
        let value = values[at];
        if (kept == 0 || values[kept - 1] != value) && elsewhere.binary_search(&value).is_err() {
            values[kept] = value;
            kept += 1;
        }
    }
    kept as u32
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
