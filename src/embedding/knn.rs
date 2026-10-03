use anyhow::Result;
use log::debug;
use ndarray::Array2;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

// Dong, Charikar and Li (2011) sample candidates at a rate on k rather than at a fixed
// count, so the cap means the same thing at every neighbour count.
const CANDIDATE_RATE: f64 = 0.5;

#[derive(Debug, Clone, Copy)]
pub struct KnnSettings {
    pub neighbours: usize,
    pub candidates: usize,
}

pub fn candidates(n_neighbours: usize) -> usize {
    ((CANDIDATE_RATE * n_neighbours as f64).ceil() as usize).max(1)
}
const MAX_ITERATIONS: usize = 20;
const MIN_WIDTH: usize = 2;

mod descent;
mod query;

pub use descent::build_knn_with;
pub use query::{nearest_exact, nearest_in};

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
pub(super) struct NeighbourList {
    dists: Vec<f64>,
    indices: Vec<u32>,
    is_new: Vec<bool>,
    fresh: Vec<bool>,
}

impl NeighbourList {
    fn new(k: usize) -> Self {
        Self {
            dists: vec![f64::INFINITY; k],
            indices: vec![u32::MAX; k],
            is_new: vec![false; k],
            fresh: vec![false; k],
        }
    }

    fn worst(&self) -> f64 {
        self.dists[self.dists.len() - 1]
    }

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
        self.fresh.copy_within(position..k - 1, position + 1);
        self.dists[position] = distance;
        self.indices[position] = index;
        self.is_new[position] = true;
        self.fresh[position] = true;
    }
}

pub(super) fn graph_of(lists: &[NeighbourList], k: usize) -> KnnGraph {
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
