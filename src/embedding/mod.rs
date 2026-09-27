pub mod features;
pub mod fuzzy;
pub mod knn;
pub mod metrics;
pub mod reach;
pub mod weight;

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

/// Ascending `indices` keep every row's columns sorted after renumbering.
pub fn lifted(graph: &Graph, indices: &[usize], rows: usize) -> Graph {
    let mut indptr = Vec::with_capacity(rows + 1);
    let mut columns = Vec::with_capacity(graph.nnz());
    let mut data = Vec::with_capacity(graph.nnz());
    indptr.push(0);
    let mut at = 0;
    for row in 0..rows {
        if indices.get(at) == Some(&row) {
            let (held, weights) = row_of(graph, at);
            columns.extend(held.iter().map(|column| indices[*column as usize] as u32));
            data.extend_from_slice(weights);
            at += 1;
        }
        indptr.push(columns.len());
    }
    sprs::CsMatI::new((rows, rows), indptr, columns, data)
}

/// The assembler's own adjacency is weak on its own, a coin flip as a must-link on the one real
/// set with a gold, so it joins the neighbour graph as another edge rather than as a constraint.
pub fn linked(
    graph: Graph,
    links: &[crate::assembly_graph::Link],
    indices: &[usize],
    weight: f32,
) -> Graph {
    let rows = graph.rows();
    let mapped = links
        .iter()
        .filter_map(|link| {
            let from = indices.binary_search(&link.from).ok()?;
            let to = indices.binary_search(&link.to).ok()?;
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
