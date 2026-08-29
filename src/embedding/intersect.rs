use sprs::TriMatI;

use crate::embedding::Graph;

const SMOOTH_TOLERANCE: f32 = 1e-5;
const SEARCH_STEPS: usize = 64;
const ABSENT_FLOOR: f32 = 1e-8;

/// Fuzzy intersection of the per-view UMAP graphs, which is what flight got from multiplying
/// its three reducers. An edge survives on the product of its strengths, so a pair has to
/// agree on every view rather than on the one view a scalar distance happens to weight.
pub fn intersect(graphs: &[Graph]) -> Graph {
    let mut result = graphs[0].clone();
    for other in &graphs[1..] {
        result = product(&result, other);
    }
    normalise_rows(&mut result);
    reset_local_metric(&mut result);
    union_with_transpose(&result)
}

fn product(left: &Graph, right: &Graph) -> Graph {
    let left_floor = absent(left);
    let right_floor = absent(right);
    let mut triplets = TriMatI::<f32, u32>::new((left.rows(), left.cols()));

    for row in 0..left.rows() {
        let (left_columns, left_data) = row_of(left, row);
        let (right_columns, right_data) = row_of(right, row);
        let mut a = 0;
        let mut b = 0;

        while a < left_columns.len() || b < right_columns.len() {
            let take_left = b >= right_columns.len()
                || (a < left_columns.len() && left_columns[a] <= right_columns[b]);
            let take_right = a >= left_columns.len()
                || (b < right_columns.len() && right_columns[b] <= left_columns[a]);

            let column = if take_left {
                left_columns[a]
            } else {
                right_columns[b]
            };
            let left_value = if take_left { left_data[a] } else { left_floor };
            let right_value = if take_right { right_data[b] } else { right_floor };

            a += usize::from(take_left);
            b += usize::from(take_right);
            triplets.add_triplet(row, column as usize, left_value * right_value);
        }
    }

    triplets.to_csr()
}

/// Half the smallest membership present rather than zero, so an edge only one view carries
/// is weakened rather than deleted.
fn absent(graph: &Graph) -> f32 {
    graph
        .data()
        .iter()
        .copied()
        .fold(f32::INFINITY, f32::min)
        .max(ABSENT_FLOOR * 2.0)
        / 2.0
}

fn row_of(graph: &Graph, row: usize) -> (&[u32], &[f32]) {
    let start = graph.indptr().index(row);
    let end = graph.indptr().index(row + 1);
    (&graph.indices()[start..end], &graph.data()[start..end])
}

fn normalise_rows(graph: &mut Graph) {
    let boundaries = (0..graph.rows())
        .map(|row| (graph.indptr().index(row), graph.indptr().index(row + 1)))
        .collect::<Vec<_>>();
    let data = graph.data_mut();
    for (start, end) in boundaries {
        let largest = data[start..end].iter().copied().fold(0.0f32, f32::max);
        if largest > 0.0 {
            for value in data[start..end].iter_mut() {
                *value /= largest;
            }
        }
    }
}

/// The product leaves each row holding less mass than the views it came from. Re-fitting the
/// local scale restores the fuzzy set cardinality UMAP's graph construction guarantees,
/// which is what `reset_local_connectivity` does before the layout sees the graph.
fn reset_local_metric(graph: &mut Graph) {
    let boundaries = (0..graph.rows())
        .map(|row| (graph.indptr().index(row), graph.indptr().index(row + 1)))
        .collect::<Vec<_>>();
    let data = graph.data_mut();

    for (start, end) in boundaries {
        let width = end - start;
        if width < 2 {
            for value in data[start..end].iter_mut() {
                *value = 1.0;
            }
            continue;
        }

        let distances = data[start..end]
            .iter()
            .map(|value| -value.max(ABSENT_FLOOR).ln())
            .collect::<Vec<_>>();
        let nearest = distances.iter().copied().fold(f32::INFINITY, f32::min);
        let scale = local_scale(&distances, nearest, (width as f32).log2());

        for (value, distance) in data[start..end].iter_mut().zip(distances.iter()) {
            *value = (-(distance - nearest).max(0.0) / scale).exp();
        }
    }
}

fn local_scale(distances: &[f32], nearest: f32, target: f32) -> f32 {
    let mut low = 0.0f32;
    let mut high = f32::INFINITY;
    let mut scale = 1.0f32;

    for _ in 0..SEARCH_STEPS {
        let total = distances
            .iter()
            .map(|distance| (-(distance - nearest).max(0.0) / scale).exp())
            .sum::<f32>();
        if (total - target).abs() < SMOOTH_TOLERANCE {
            break;
        }
        if total > target {
            high = scale;
            scale = (low + high) / 2.0;
        } else {
            low = scale;
            scale = if high.is_infinite() {
                scale * 2.0
            } else {
                (low + high) / 2.0
            };
        }
    }

    scale.max(ABSENT_FLOOR)
}

fn union_with_transpose(graph: &Graph) -> Graph {
    let transposed = graph.transpose_view().to_csr();
    let mut triplets = TriMatI::<f32, u32>::new((graph.rows(), graph.cols()));

    for row in 0..graph.rows() {
        let (columns, data) = row_of(graph, row);
        let (other_columns, other_data) = row_of(&transposed, row);
        let mut a = 0;
        let mut b = 0;

        while a < columns.len() || b < other_columns.len() {
            let take_left =
                b >= other_columns.len() || (a < columns.len() && columns[a] <= other_columns[b]);
            let take_right =
                a >= columns.len() || (b < other_columns.len() && other_columns[b] <= columns[a]);

            let column = if take_left { columns[a] } else { other_columns[b] };
            let left = if take_left { data[a] } else { 0.0 };
            let right = if take_right { other_data[b] } else { 0.0 };

            a += usize::from(take_left);
            b += usize::from(take_right);
            let value = left + right - left * right;
            if value > 0.0 {
                triplets.add_triplet(row, column as usize, value);
            }
        }
    }

    triplets.to_csr()
}
