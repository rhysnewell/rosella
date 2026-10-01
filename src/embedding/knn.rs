use anyhow::Result;
use log::debug;
use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
    sync::{
        Mutex,
        atomic::{AtomicU64, Ordering},
    },
};

// Dong, Charikar and Li (2011) sample candidates at a rate on k rather than at a fixed
// count, so the cap means the same thing at every neighbour count.
const CANDIDATE_RATE: f64 = 0.5;

pub fn candidates(n_neighbours: usize) -> usize {
    ((CANDIDATE_RATE * n_neighbours as f64).ceil() as usize).max(1)
}
const MAX_ITERATIONS: usize = 20;
const MIN_WIDTH: usize = 2;

pub trait Metric: Sync {
    fn distance(&self, a: usize, b: usize) -> f64;

    // Most pairs the descent meets cannot enter either list, and a metric that can show that
    // cheaply returns any value past `bound` instead of the distance.
    fn within(&self, a: usize, b: usize, _bound: f64) -> f64 {
        self.distance(a, b)
    }

    // Pairs sharing a side are measured together so a metric can overlap their arithmetic.
    fn within_many(
        &self,
        a: usize,
        others: &[u32],
        bound: impl Fn(usize) -> f64,
        out: &mut Vec<f64>,
    ) {
        out.extend(
            others
                .iter()
                .map(|b| self.within(a, *b as usize, bound(*b as usize))),
        );
    }
}

impl<F: Fn(usize, usize) -> f64 + Sync> Metric for F {
    fn distance(&self, a: usize, b: usize) -> f64 {
        self(a, b)
    }
}

pub struct KnnGraph {
    pub indices: Array2<u32>,
    pub dists: Array2<f32>,
}

impl KnnGraph {
    pub fn n_points(&self) -> usize {
        self.indices.nrows()
    }

    // Columns are sorted by distance then index, so the first k of a wider build are exactly
    // the k nearest and a sparser round needs no build of its own.
    pub fn truncate(&self, k: usize) -> KnnGraph {
        let width = k.min(self.indices.ncols());
        KnnGraph {
            indices: self.indices.slice(ndarray::s![.., ..width]).to_owned(),
            dists: self.dists.slice(ndarray::s![.., ..width]).to_owned(),
        }
    }

    // A surviving row is exact for any width up to its own survivor count, so the width is the
    // narrowest row and a shorter one cannot be padded past the manifold builders.
    pub fn induced(&self, keep: &[usize]) -> Option<KnnGraph> {
        let survivors = keep
            .iter()
            .map(|node| {
                self.indices
                    .row(*node)
                    .iter()
                    .filter(|column| keep.binary_search(&(**column as usize)).is_ok())
                    .count()
            })
            .collect::<Vec<_>>();
        let width = survivors.iter().copied().min()?;
        if width < MIN_WIDTH {
            if log::log_enabled!(log::Level::Debug) {
                let mut spread = survivors;
                spread.sort_unstable();
                let thin = spread
                    .iter()
                    .take_while(|count| **count < MIN_WIDTH)
                    .count();
                debug!(
                    "Induced {} rows refused: {thin} under {MIN_WIDTH}, tenth {}, median {}",
                    keep.len(),
                    spread[spread.len() / 10],
                    spread[spread.len() / 2]
                );
            }
            return None;
        }
        let mut indices = Array2::<u32>::zeros((keep.len(), width));
        let mut dists = Array2::<f32>::zeros((keep.len(), width));
        for (row, node) in keep.iter().enumerate() {
            let mut taken = 0;
            for (column, distance) in self.indices.row(*node).iter().zip(self.dists.row(*node)) {
                if taken == width {
                    break;
                }
                let Ok(position) = keep.binary_search(&(*column as usize)) else {
                    continue;
                };
                indices[[row, taken]] = position as u32;
                dists[[row, taken]] = *distance;
                taken += 1;
            }
        }
        Some(KnnGraph { indices, dists })
    }
}

pub fn write_report(views: &[(&str, KnnGraph)], names: &[&str], path: &Path) -> Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "view\tcontig\trank\tneighbour\tdistance")?;
    for (view, knn) in views {
        for row in 0..knn.n_points() {
            for (rank, neighbour) in knn.indices.row(row).iter().enumerate() {
                if *neighbour == u32::MAX {
                    continue;
                }
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    view,
                    names[row],
                    rank,
                    names[*neighbour as usize],
                    knn.dists[[row, rank]]
                )?;
            }
        }
    }
    out.flush()?;
    Ok(())
}

// Ties break on index, so insertion is order independent and the build reproducible under rayon.
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

    fn worst(&self) -> f64 {
        self.dists[self.dists.len() - 1]
    }

    // Total order on (distance, index). The index tiebreak is what keeps insertion
    // order independent.
    fn sorts_before(&self, position: usize, distance: f64, index: u32) -> bool {
        self.dists[position] < distance
            || (self.dists[position] == distance && self.indices[position] < index)
    }

    fn push(&mut self, distance: f64, index: u32) {
        let k = self.dists.len();
        if self.sorts_before(k - 1, distance, index) || self.indices.contains(&index) {
            return;
        }

        let (mut position, mut above) = (0, k - 1);
        while position < above {
            let middle = (position + above) / 2;
            match self.sorts_before(middle, distance, index) {
                true => position = middle + 1,
                false => above = middle,
            }
        }
        self.dists.copy_within(position..k - 1, position + 1);
        self.indices.copy_within(position..k - 1, position + 1);
        self.is_new.copy_within(position..k - 1, position + 1);
        self.dists[position] = distance;
        self.indices[position] = index;
        self.is_new[position] = true;
    }
}

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
    descend(lists, n, k, max_candidates, metric)
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
    metric: M,
) -> KnnGraph {
    let progress = crate::progress::counted(
        crate::progress::Stage::NearestNeighbours,
        MAX_ITERATIONS as u64,
    );
    let mut candidates = Candidates::new(n, max_candidates.max(1));
    for round in 0..MAX_ITERATIONS {
        build_candidates(&lists.rows, n, &mut candidates);

        (0..n)
            .into_par_iter()
            .for_each(|i| join(&metric, &lists, candidates.new_of(i), candidates.old_of(i)));

        let taken = taken_slots(&lists.rows);
        debug!("Descent round {round} took {taken} slots");
        progress.inc(1);
        progress.set_message(format!("{taken} neighbours taken"));
        if taken as f64 <= crate::tuning::CONVERGENCE_FRACTION * k as f64 * n as f64 {
            break;
        }
    }
    progress.finish_and_clear();

    let mut indices = Array2::from_elem((n, k), u32::MAX);
    let mut dists = Array2::from_elem((n, k), f32::INFINITY);
    for (i, list) in lists.rows.iter().enumerate() {
        let list = list.lock().unwrap();
        for j in 0..k {
            indices[[i, j]] = list.indices[j];
            dists[[i, j]] = list.dists[j] as f32;
        }
    }
    KnnGraph { indices, dists }
}

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
pub fn nearest_exact<M: Metric>(base: usize, queries: &[usize], k: usize, metric: M) -> KnnGraph {
    let k = k.min(base.saturating_sub(1)).max(1);
    let lists = queries
        .par_iter()
        .map(|query| {
            let mut list = NeighbourList::new(k);
            let (mut chunk, mut distances) = (Vec::with_capacity(k), Vec::with_capacity(k));
            for start in (0..base).step_by(k) {
                chunk.clear();
                distances.clear();
                chunk.extend(
                    (start..(start + k).min(base))
                        .filter(|other| other != query)
                        .map(|other| other as u32),
                );
                let bound = list.worst();
                metric.within_many(*query, &chunk, |_| bound, &mut distances);
                for (other, distance) in chunk.iter().zip(&distances) {
                    list.push(*distance, *other);
                }
            }
            list
        })
        .collect::<Vec<_>>();
    graph_of(&lists, k)
}

fn graph_of(lists: &[NeighbourList], k: usize) -> KnnGraph {
    let mut indices = Array2::from_elem((lists.len(), k), u32::MAX);
    let mut dists = Array2::from_elem((lists.len(), k), f32::INFINITY);
    for (row, list) in lists.iter().enumerate() {
        for slot in 0..k {
            indices[[row, slot]] = list.indices[slot];
            dists[[row, slot]] = list.dists[slot] as f32;
        }
    }
    KnnGraph { indices, dists }
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

// Forward and reverse neighbour lists, split by whether the edge is new since the last
// pass. Every bucket is capped, so one flat allocation of that stride holds the whole pass.
struct Candidates {
    stride: usize,
    new: Vec<u32>,
    new_len: Vec<u32>,
    old: Vec<u32>,
    old_len: Vec<u32>,
}

fn push_candidate(
    values: &mut [u32],
    lengths: &mut [u32],
    stride: usize,
    row: usize,
    value: u32,
) -> bool {
    let length = lengths[row] as usize;
    if length == stride {
        return false;
    }
    values[row * stride + length] = value;
    lengths[row] = length as u32 + 1;
    true
}

impl Candidates {
    fn new(n: usize, stride: usize) -> Self {
        Self {
            stride,
            new: vec![u32::MAX; n * stride],
            new_len: vec![0; n],
            old: vec![u32::MAX; n * stride],
            old_len: vec![0; n],
        }
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

// Built in index order so the caps fall the same way every run.
fn build_candidates(neighbours: &[Mutex<NeighbourList>], n: usize, candidates: &mut Candidates) {
    candidates.clear();
    let stride = candidates.stride;

    for (i, entry) in neighbours.iter().enumerate().take(n) {
        let mut list = entry.lock().unwrap();
        for slot in 0..list.indices.len() {
            let neighbour = list.indices[slot];
            if neighbour == u32::MAX {
                continue;
            }
            let (values, lengths) = if list.is_new[slot] {
                (&mut candidates.new, &mut candidates.new_len)
            } else {
                (&mut candidates.old, &mut candidates.old_len)
            };
            if push_candidate(values, lengths, stride, i, neighbour) {
                list.is_new[slot] = false;
            }
            push_candidate(values, lengths, stride, neighbour as usize, i as u32);
        }
    }
    candidates.dedup();
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
            lists.offer(a, *distance, *b);
            lists.offer(*b as usize, *distance, a as u32);
        }
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
                .is_new
                .iter()
                .filter(|fresh| **fresh)
                .count()
        })
        .sum()
}
