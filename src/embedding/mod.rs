pub mod features;
pub mod knn;
pub mod metrics;
pub mod fuzzy;

pub const KNN_ASSEMBLY: &str = "knn_assembly";
pub const KNN_POOL: &str = "knn_pool";
pub const KNN_SPLIT: &str = "knn_split";

pub type Graph = sprs::CsMatI<f32, u32, usize>;

pub(crate) fn row_of(graph: &Graph, row: usize) -> (&[u32], &[f32]) {
    let start = graph.indptr().index(row);
    let end = graph.indptr().index(row + 1);
    (&graph.indices()[start..end], &graph.data()[start..end])
}

/// `Level::from_graph` and `label_propagation` size their vectors by `graph.rows()` and index
/// them by node id, so a slice has to come back renumbered rather than masked.
pub fn induced(graph: &Graph, indices: &[usize]) -> Graph {
    let mut triplets = sprs::TriMatI::<f32, u32>::new((indices.len(), indices.len()));
    for (row, node) in indices.iter().enumerate() {
        let (columns, weights) = row_of(graph, *node);
        for (column, weight) in columns.iter().zip(weights) {
            if let Ok(position) = indices.binary_search(&(*column as usize)) {
                triplets.add_triplet(row, position, *weight);
            }
        }
    }
    triplets.to_csr()
}

/// The assembler's own adjacency is weak on its own, a coin flip as a must-link on the one real
/// set with a gold, so it joins the neighbour graph as another edge rather than as a constraint.
pub fn linked(graph: Graph, links: &[(usize, usize)], indices: &[usize], weight: f32) -> Graph {
    let rows = graph.rows();
    let mapped = links
        .iter()
        .filter_map(|(from, to)| {
            let from = indices.binary_search(from).ok()?;
            let to = indices.binary_search(to).ok()?;
            (from != to).then_some((from.min(to), from.max(to)))
        })
        .collect::<std::collections::HashSet<_>>();
    if mapped.is_empty() {
        return graph;
    }
    let mut triplets = sprs::TriMatI::<f32, u32>::new((rows, rows));
    let mut held = std::collections::HashSet::new();
    for row in 0..rows {
        let (columns, weights) = row_of(&graph, row);
        for (column, edge) in columns.iter().zip(weights) {
            let column = *column as usize;
            let pair = (row.min(column), row.max(column));
            let raised = match mapped.contains(&pair) {
                true => edge.max(weight),
                false => *edge,
            };
            triplets.add_triplet(row, column, raised);
            held.insert(pair);
        }
    }
    for (from, to) in mapped.difference(&held) {
        triplets.add_triplet(*from, *to, weight);
        triplets.add_triplet(*to, *from, weight);
    }
    triplets.to_csr()
}
