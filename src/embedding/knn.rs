use anyhow::Result;
use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
    sync::Mutex,
};

pub const MAX_CANDIDATES: usize = 32;
const MAX_ITERATIONS: usize = 20;
const CONVERGENCE_FRACTION: f64 = 0.001;
const ROW_SEED_STRIDE: u64 = 0x9E37_79B9_7F4A_7C15;
const MIN_WIDTH: usize = 2;


pub struct KnnGraph {
    pub indices: Array2<u32>,
    pub dists: Array2<f32>,
}

impl KnnGraph {
    pub fn n_points(&self) -> usize {
        self.indices.nrows()
    }

    /// Columns are sorted by distance then index, so the first k of a wider build are exactly
    /// the k nearest and a sparser round needs no build of its own.
    pub fn truncate(&self, k: usize) -> KnnGraph {
        let width = k.min(self.indices.ncols());
        KnnGraph {
            indices: self.indices.slice(ndarray::s![.., ..width]).to_owned(),
            dists: self.dists.slice(ndarray::s![.., ..width]).to_owned(),
        }
    }

    /// A surviving row is exact for any width up to its own survivor count, so the width is the
    /// narrowest row and a shorter one cannot be padded past the manifold builders.
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
        let width = survivors.iter().copied().min().unwrap_or(0);
        if width < MIN_WIDTH {
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

/// Deterministic for a given seed and k, whatever the thread count, because the only
/// randomness is the seeded initial sample and every later step is order independent.
pub fn build_knn_with<M>(
    n: usize,
    k: usize,
    max_candidates: usize,
    seed: u64,
    metric: M,
) -> KnnGraph
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

    let mut candidates = Candidates::new(n, max_candidates.max(1));
    for _ in 0..MAX_ITERATIONS {
        build_candidates(&neighbours, n, &mut candidates);

        let updates: usize = (0..n)
            .into_par_iter()
            .map(|i| {
                join(
                    &metric,
                    &neighbours,
                    candidates.new_of(i),
                    candidates.old_of(i),
                )
            })
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
/// pass. Every bucket is capped, so one flat allocation of that stride holds the whole pass.
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
}

/// Built in index order so the caps fall the same way every run.
fn build_candidates(neighbours: &[Mutex<NeighbourList>], n: usize, candidates: &mut Candidates) {
    candidates.clear();
    let stride = candidates.stride;

    for i in 0..n {
        let mut list = neighbours[i].lock().unwrap();
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
