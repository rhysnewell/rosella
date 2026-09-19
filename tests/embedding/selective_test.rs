use ndarray::Array2;
use rosella::embedding::fuzzy::{manifold_graph, manifold_graph_with, scales};
use rosella::embedding::knn::KnnGraph;
use rosella::embedding::selective::{hop_extras, needy, wide_extras};

fn graph_of(rows: &[(Vec<u32>, Vec<f32>)]) -> KnnGraph {
    let width = rows[0].0.len();
    let mut indices = Array2::<u32>::zeros((rows.len(), width));
    let mut dists = Array2::<f32>::zeros((rows.len(), width));
    for (row, (neighbours, distances)) in rows.iter().enumerate() {
        for column in 0..width {
            indices[[row, column]] = neighbours[column];
            dists[[row, column]] = distances[column];
        }
    }
    KnnGraph { indices, dists }
}

fn ring(points: usize, width: usize) -> KnnGraph {
    let rows = (0..points)
        .map(|point| {
            let neighbours = (1..=width)
                .map(|step| ((point + step) % points) as u32)
                .collect::<Vec<_>>();
            let distances = (1..=width).map(|step| step as f32 * 0.1).collect();
            (neighbours, distances)
        })
        .collect::<Vec<_>>();
    graph_of(&rows)
}

/// The whole premise of the selective widen: a build carrying spare columns has to leave every
/// contig nobody selected exactly where a build of the base width would have left it.
#[test]
fn spare_columns_change_nothing_until_a_contig_is_selected() {
    let wide = ring(40, 12);
    let base = 4;
    let narrow = wide.truncate(base);

    let held = manifold_graph(40, &wide, base);
    let shipped = manifold_graph(40, &narrow, base);
    assert_eq!(held.indices(), shipped.indices());
    assert_eq!(held.data(), shipped.data());

    let (sigmas, rhos) = scales(wide.dists.view(), base);
    let mut chosen = vec![false; 40];
    chosen[7] = true;
    let extras = wide_extras(&wide, base, &chosen);
    let widened = manifold_graph_with(40, &wide, base, &sigmas, &rhos, &extras);
    assert!(widened.nnz() > shipped.nnz());
}

#[test]
fn a_selected_row_carries_the_columns_past_the_base_width_and_no_other_row_does() {
    let wide = ring(10, 6);
    let extras = wide_extras(
        &wide,
        2,
        &[
            true, false, false, false, false, false, false, false, false, false,
        ],
    );
    assert_eq!(extras[0], vec![(3, 0.3), (4, 0.4), (5, 0.5), (6, 0.6)]);
    assert!(extras[1..].iter().all(Vec::is_empty));
}

#[test]
fn the_second_hop_offers_only_what_the_first_hop_missed() {
    let knn = graph_of(&[
        (vec![1, 2], vec![0.1, 0.2]),
        (vec![0, 3], vec![0.1, 0.4]),
        (vec![0, 1], vec![0.2, 0.3]),
        (vec![1, 2], vec![0.4, 0.3]),
    ]);
    let extras = hop_extras(&knn, 2, 5, &[true, false, false, false], |a, b| {
        (a as f64 - b as f64).abs()
    });
    assert_eq!(extras[0], vec![(3, 3.0)]);
    assert!(extras[1..].iter().all(Vec::is_empty));
}

#[test]
fn the_quantile_selects_the_share_and_keeps_a_tie_whole() {
    assert_eq!(
        needy(&[0.9, 0.1, 0.5, 0.2], 0.25),
        vec![true, false, false, false]
    );
    assert_eq!(
        needy(&[0.9, 0.9, 0.1, 0.1], 0.25),
        vec![true, true, false, false]
    );
    assert_eq!(needy(&[0.9, 0.1, 0.5, 0.2], 0.0), vec![false; 4]);
}
