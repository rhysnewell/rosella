use ndarray::ArrayView2;
use rayon::prelude::*;
use sprs::CsMatI;

use crate::embedding::Graph;
use crate::embedding::knn::KnnGraph;

// McInnes, Healy & Melville (2018), section 3.1, for the graph alone and not the layout.

const TOLERANCE: f32 = 1e-5;
const MIN_SCALE: f32 = 1e-3;
const ITERATIONS: usize = 64;

// At 1.0 every contig keeps an edge to its nearest neighbour whatever the distance, which stops a
// sparse region fragmenting. 2.0 was measured over 14 sets and lost t1 on nine.
const LOCAL_CONNECTIVITY: f32 = 1.0;

// A pure fuzzy union keeps an edge either direction found. Discounting one-sided edges at 0.5 lost,
// since a fragment its neighbour does not reciprocate still belongs in that genome.
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

// Binary search for the bandwidth that makes each contig's memberships sum to log2(k), so a
// contig in a dense region and one in a sparse region contribute comparable edges.
fn bandwidth(row: &[f32], rho: f32, target: f32, floor: f32) -> f32 {
    let (mut low, mut high, mut mid) = (0.0f32, f32::INFINITY, 1.0f32);
    for _ in 0..ITERATIONS {
        let mut total = 0.0;
        // The rows hold no self entry, so this leaves the nearest neighbour out and widens every
        // sigma. Summing from it was measured and lost bins on every family.
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

// A row the search could not fill is padded with infinity, which would make the floor and so
// every membership in the row infinite and flat.
fn finite_mean<'a>(distances: impl Iterator<Item = &'a f32>) -> f32 {
    let (total, count) = distances
        .filter(|distance| distance.is_finite())
        .fold((0.0f32, 0usize), |(total, count), distance| {
            (total + distance, count + 1)
        });
    if count == 0 {
        0.0
    } else {
        total / count as f32
    }
}

pub fn scales(distances: ArrayView2<f32>, k: usize) -> (Vec<f32>, Vec<f32>) {
    let width = k.min(distances.ncols());
    let target = (width as f32).log2();
    (0..distances.nrows())
        .into_par_iter()
        .map(|point| {
            let row = distances.row(point);
            let row = &row.as_slice().expect("knn distances are contiguous")[..width];
            let rho = nearest(row, LOCAL_CONNECTIVITY);
            let floor = MIN_SCALE * finite_mean(row.iter());
            (bandwidth(row, rho, target, floor), rho)
        })
        .unzip()
}

// Rows are independent, so each is built whole and the pieces concatenated. That keeps the
// index arithmetic out of the parallel section.
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

// An edge is directed until here: contig A may hold B as a neighbour without B holding A.
// The union keeps either direction, which is what makes the graph symmetric.
fn union(graph: &Graph) -> Graph {
    let points = graph.rows();
    let product = 1.0 - 2.0 * SET_OP_MIX;

    let mut incoming = vec![Vec::new(); points];
    for row in 0..points {
        let (start, end) = (graph.indptr().index(row), graph.indptr().index(row + 1));
        for (column, weight) in graph.indices()[start..end]
            .iter()
            .zip(&graph.data()[start..end])
        {
            incoming[*column as usize].push((row as u32, *weight));
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
            for (source, forward) in &incoming[point] {
                if at(graph, point, *source) != 0.0 {
                    continue;
                }
                let value = SET_OP_MIX * forward;
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

pub fn manifold_graph(points: usize, knn: &KnnGraph) -> Graph {
    let _timer = crate::timing::scope("manifold");
    let width = knn.indices.ncols();
    let (sigmas, rhos) = scales(knn.dists.view(), width);
    union(&memberships(points, knn, width, &sigmas, &rhos))
}
