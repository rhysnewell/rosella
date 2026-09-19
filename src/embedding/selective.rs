use std::cmp::Ordering;
use std::collections::HashMap;

use ndarray::ArrayView2;
use rayon::prelude::*;

use crate::embedding::fuzzy::membership;
use crate::embedding::knn::KnnGraph;

/// rsl-split-graph-reach put the nearest same-genome contig of an unreachable contig at a median
/// rank of 208, so three times the shipped hundred covers the median without another knob.
const WIDE_FACTOR: usize = 3;

pub type Extras = Vec<Vec<(u32, f32)>>;

#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum Mode {
    #[default]
    Off,
    Wide,
    Hop,
}

impl Mode {
    pub fn parse(name: &str) -> Option<Mode> {
        match name.trim() {
            "off" => Some(Mode::Off),
            "wide" => Some(Mode::Wide),
            "hop" => Some(Mode::Hop),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Reach {
    pub mode: Mode,
    pub share: f64,
}

impl Reach {
    pub fn build_width(&self, base: usize) -> usize {
        match self.mode {
            Mode::Wide => base * WIDE_FACTOR,
            _ => base,
        }
    }
}

pub fn budget(base: usize) -> usize {
    base * (WIDE_FACTOR - 1)
}

/// A neighbourhood is exhausted when its own outermost edge still carries a high membership: the
/// hundredth neighbour is barely further than the first, so the cut at a hundred is arbitrary.
pub fn edge_membership(
    dists: ArrayView2<f32>,
    width: usize,
    sigmas: &[f32],
    rhos: &[f32],
) -> Vec<f32> {
    let last = width.min(dists.ncols()).saturating_sub(1);
    (0..dists.nrows())
        .into_par_iter()
        .map(|point| membership(dists[(point, last)], rhos[point], sigmas[point]))
        .collect()
}

/// A quantile rather than a bar, so the rule tracks the assembly's own spread instead of a
/// constant measured on one benchmark.
pub fn needy(edge: &[f32], share: f64) -> Vec<bool> {
    let take = (edge.len() as f64 * share).round() as usize;
    if take == 0 || edge.is_empty() {
        return vec![false; edge.len()];
    }
    let mut sorted = edge.to_vec();
    sorted.sort_unstable_by(|a, b| b.partial_cmp(a).unwrap_or(Ordering::Equal));
    let bar = sorted[take.min(sorted.len()) - 1];
    edge.iter().map(|value| *value >= bar).collect()
}

pub fn wide_extras(knn: &KnnGraph, width: usize, chosen: &[bool]) -> Extras {
    let columns = knn.indices.ncols();
    (0..knn.n_points())
        .into_par_iter()
        .map(|point| match chosen.get(point) {
            Some(true) => (width..columns)
                .map(|position| (knn.indices[(point, position)], knn.dists[(point, position)]))
                .collect(),
            _ => Vec::new(),
        })
        .collect()
}

/// The two hop sum bounds the real distance under the triangle inequality, so it shortlists for
/// free and only the shortlist pays for a metric call.
pub fn hop_extras(
    knn: &KnnGraph,
    width: usize,
    budget: usize,
    chosen: &[bool],
    distance: impl Fn(usize, usize) -> f64 + Sync,
) -> Extras {
    let width = width.min(knn.indices.ncols());
    (0..knn.n_points())
        .into_par_iter()
        .map(|point| {
            if chosen.get(point) != Some(&true) {
                return Vec::new();
            }
            let own = knn.indices.row(point);
            let held = own.iter().take(width).copied().collect::<Vec<_>>();
            let mut bounds: HashMap<u32, f32> = HashMap::new();
            for position in 0..width {
                let near = knn.indices[(point, position)] as usize;
                if near == point || near >= knn.n_points() {
                    continue;
                }
                let first = knn.dists[(point, position)];
                for step in 0..width {
                    let far = knn.indices[(near, step)];
                    if far as usize == point || held.contains(&far) {
                        continue;
                    }
                    let bound = first + knn.dists[(near, step)];
                    bounds
                        .entry(far)
                        .and_modify(|held| *held = held.min(bound))
                        .or_insert(bound);
                }
            }
            let mut shortlist = bounds.into_iter().collect::<Vec<_>>();
            shortlist.sort_unstable_by(|a, b| {
                a.1.partial_cmp(&b.1)
                    .unwrap_or(Ordering::Equal)
                    .then(a.0.cmp(&b.0))
            });
            shortlist.truncate(budget);
            let mut scored = shortlist
                .into_iter()
                .map(|(far, _)| (far, distance(point, far as usize) as f32))
                .collect::<Vec<_>>();
            scored.sort_unstable_by(|a, b| {
                a.1.partial_cmp(&b.1)
                    .unwrap_or(Ordering::Equal)
                    .then(a.0.cmp(&b.0))
            });
            scored
        })
        .collect()
}
