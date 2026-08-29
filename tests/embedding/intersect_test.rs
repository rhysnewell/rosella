use rosella::embedding::Graph;
use rosella::embedding::intersect::intersect;
use sprs::TriMatI;

const SHARED: [(usize, usize, f32); 13] = [
    (0, 1, 0.9),
    (0, 2, 0.8),
    (0, 4, 0.3),
    (0, 5, 0.25),
    (1, 2, 0.7),
    (1, 3, 0.45),
    (1, 4, 0.35),
    (1, 5, 0.2),
    (2, 3, 0.5),
    (2, 4, 0.4),
    (2, 5, 0.3),
    (3, 4, 0.95),
    (3, 5, 0.55),
];

fn view(pair: Option<f32>) -> Graph {
    let mut triplets = TriMatI::<f32, u32>::new((6, 6));
    let mut edges = SHARED.to_vec();
    if let Some(strength) = pair {
        edges.push((0, 3, strength));
    }
    for (row, column, value) in edges {
        triplets.add_triplet(row, column, value);
        triplets.add_triplet(column, row, value);
    }
    triplets.to_csr()
}

fn weight(graph: &Graph, row: usize, column: usize) -> f32 {
    graph.get(row, column).copied().unwrap_or(0.0)
}

/// What a scalar distance cannot express. The pair 0,3 is close on coverage either way; the
/// second view's opinion is what decides how close they end up. A weighted product of
/// distances lets the coverage term carry the pair on its own.
#[test]
fn the_second_view_moves_an_edge_the_first_agrees_on() {
    let coverage = view(Some(0.85));
    let agrees = intersect(&[coverage.clone(), view(Some(0.8))]);
    let disagrees = intersect(&[coverage, view(Some(0.1))]);

    assert!(
        weight(&agrees, 0, 3) > weight(&disagrees, 0, 3),
        "agreement scored {} against {} for disagreement",
        weight(&agrees, 0, 3),
        weight(&disagrees, 0, 3)
    );
}

/// An edge one view never saw is weakened, not deleted, so a contig whose views disagree
/// everywhere still reaches the graph rather than becoming a point the layout cannot place.
#[test]
fn an_edge_in_one_view_alone_survives() {
    let result = intersect(&[view(Some(0.85)), view(None)]);

    assert!(weight(&result, 0, 3) > 0.0);
    assert!(weight(&result, 0, 3) < weight(&result, 0, 1));
}

/// The layout reads the graph as undirected, and the per-row rescaling before the union is
/// what would otherwise leave it lopsided.
#[test]
fn the_result_is_symmetric() {
    let result = intersect(&[view(Some(0.85)), view(Some(0.2))]);

    for row in 0..6 {
        for column in 0..6 {
            let forward = weight(&result, row, column);
            let backward = weight(&result, column, row);
            assert!(
                (forward - backward).abs() < 1e-6,
                "{row},{column} held {forward} against {backward}"
            );
        }
    }
}
