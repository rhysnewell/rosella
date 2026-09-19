use ndarray::ArrayView2;
use rayon::prelude::*;
use sprs::CsMatI;

use crate::embedding::Graph;
use crate::embedding::knn::KnnGraph;

/// McInnes, Healy & Melville (2018), section 3.1. Only the graph is used here; the low
/// dimensional layout that method is better known for is not part of this pipeline.
const TOLERANCE: f32 = 1e-5;
const MIN_SCALE: f32 = 1e-3;
const ITERATIONS: usize = 64;

/// Neighbours assumed connected at a local level. At 1.0 every contig keeps an edge to its
/// nearest neighbour whatever the distance, which is what stops a sparse region fragmenting.
/// 2.0 was measured over 14 sets and lost t1 on nine of them.
const LOCAL_CONNECTIVITY: f32 = 1.0;

/// A pure fuzzy union, so an edge either direction found survives whole. Discounting the
/// one-sided edges at 0.5 was measured and lost, so a fragmented contig its neighbour does
/// not reciprocate still belongs in that genome.
const SET_OP_MIX: f32 = 1.0;

fn nearest(row: &[f32], connectivity: f32) -> f32 {
    let index = connectivity.floor() as usize;
    let interpolation = connectivity - connectivity.floor();
    let non_zero = row.iter().copied().filter(|distance| *distance > 0.0);

    if row.iter().filter(|distance| **distance > 0.0).count() < index.max(1) {
        return row
            .iter()
            .copied()
            .filter(|distance| *distance > 0.0)
            .fold(0.0f32, f32::max);
    }
    if index == 0 {
        return interpolation * non_zero.clone().next().unwrap_or(0.0);
    }
    let mut held = non_zero.skip(index - 1);
    let previous = held.next().unwrap_or(0.0);
    match interpolation > TOLERANCE {
        true => previous + interpolation * (held.next().unwrap_or(0.0) - previous),
        false => previous,
    }
}

/// Binary search for the bandwidth that makes each contig's memberships sum to log2(k), so a
/// contig in a dense region and one in a sparse region contribute comparable edges.
fn bandwidth(row: &[f32], rho: f32, target: f32, floor: f32) -> f32 {
    let (mut low, mut high, mut mid) = (0.0f32, f32::INFINITY, 1.0f32);
    for _ in 0..ITERATIONS {
        let mut total = 0.0;
        for distance in row.iter().skip(1) {
            let gap = distance - rho;
            total += match gap > 0.0 {
                true => (-(gap / mid)).exp(),
                false => 1.0,
            };
        }
        if (total - target).abs() < TOLERANCE {
            break;
        }
        if total > target {
            high = mid;
            mid = (low + high) / 2.0;
        } else {
            low = mid;
            mid = match high.is_infinite() {
                true => mid * 2.0,
                false => (low + high) / 2.0,
            };
        }
    }
    mid.max(floor)
}

pub(crate) fn membership(distance: f32, rho: f32, sigma: f32) -> f32 {
    match distance - rho <= 0.0 || sigma == 0.0 {
        true => 1.0,
        false => (-((distance - rho) / sigma)).exp(),
    }
}

pub fn scales(distances: ArrayView2<f32>, k: usize) -> (Vec<f32>, Vec<f32>) {
    let width = k.min(distances.ncols());
    let target = (width as f32).log2();
    let held = distances.slice(ndarray::s![.., ..width]);
    let overall = held.mean().unwrap_or(0.0);
    (0..distances.nrows())
        .into_par_iter()
        .map(|point| {
            let row = distances.row(point);
            let row = &row.as_slice().expect("knn distances are contiguous")[..width];
            let rho = nearest(row, LOCAL_CONNECTIVITY);
            let mean = row.iter().sum::<f32>() / row.len() as f32;
            let floor = MIN_SCALE * if rho > 0.0 { mean } else { overall };
            (bandwidth(row, rho, target, floor), rho)
        })
        .unzip()
}

/// Rows are independent, so each is built whole and the pieces concatenated. That keeps the
/// index arithmetic out of the parallel section.
fn rows_into_graph(points: usize, rows: Vec<Vec<(u32, f32)>>) -> Graph {
    let mut indptr = Vec::with_capacity(points + 1);
    let mut indices = Vec::new();
    let mut data = Vec::new();
    indptr.push(0);
    for row in rows {
        for (column, value) in row {
            indices.push(column);
            data.push(value);
        }
        indptr.push(indices.len());
    }
    CsMatI::new((points, points), indptr, indices, data)
}

fn memberships(points: usize, knn: &KnnGraph, width: usize, sigmas: &[f32], rhos: &[f32]) -> Graph {
    let rows = (0..points)
        .into_par_iter()
        .map(|point| {
            let mut row = (0..width)
                .map(|position| (knn.indices[(point, position)], knn.dists[(point, position)]))
                .filter_map(|(neighbour, distance)| {
                    if neighbour as usize == point || neighbour as usize >= points {
                        return None;
                    }
                    let value = membership(distance, rhos[point], sigmas[point]);
                    (value != 0.0).then_some((neighbour, value))
                })
                .collect::<Vec<_>>();
            row.sort_unstable_by_key(|(column, _)| *column);
            row.dedup_by_key(|(column, _)| *column);
            row
        })
        .collect::<Vec<_>>();
    rows_into_graph(points, rows)
}

fn at(graph: &Graph, row: usize, column: u32) -> f32 {
    let (start, end) = (graph.indptr().index(row), graph.indptr().index(row + 1));
    match graph.indices()[start..end].binary_search(&column) {
        Ok(found) => graph.data()[start + found],
        Err(_) => 0.0,
    }
}

/// An edge is directed until here: contig A may hold B as a neighbour without B holding A.
/// The union keeps either direction, which is what makes the graph symmetric.
fn union(graph: &Graph) -> Graph {
    let points = graph.rows();
    let product = 1.0 - 2.0 * SET_OP_MIX;

    let mut incoming = vec![Vec::new(); points];
    for row in 0..points {
        let (start, end) = (graph.indptr().index(row), graph.indptr().index(row + 1));
        for column in &graph.indices()[start..end] {
            incoming[*column as usize].push(row as u32);
        }
    }

    let rows = (0..points)
        .into_par_iter()
        .map(|point| {
            let (start, end) = (graph.indptr().index(point), graph.indptr().index(point + 1));
            let mut row = Vec::with_capacity(end - start + incoming[point].len());
            for (column, forward) in graph.indices()[start..end]
                .iter()
                .zip(&graph.data()[start..end])
            {
                let back = at(graph, *column as usize, point as u32);
                let value = SET_OP_MIX * forward + SET_OP_MIX * back + product * forward * back;
                if value != 0.0 {
                    row.push((*column, value));
                }
            }
            for source in &incoming[point] {
                if at(graph, point, *source) != 0.0 {
                    continue;
                }
                let value = SET_OP_MIX * at(graph, *source as usize, point as u32);
                if value != 0.0 {
                    row.push((*source, value));
                }
            }
            row.sort_unstable_by_key(|(column, _)| *column);
            row
        })
        .collect::<Vec<_>>();
    rows_into_graph(points, rows)
}

pub fn manifold_graph(points: usize, knn: &KnnGraph, n_neighbours: usize) -> Graph {
    let _timer = crate::timing::scope("manifold");
    let width = n_neighbours.min(knn.indices.ncols());
    let (sigmas, rhos) = scales(knn.dists.view(), width);
    union(&memberships(points, knn, width, &sigmas, &rhos))
}
