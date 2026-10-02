pub mod features;
pub mod fuzzy;
pub mod knn;
pub mod metrics;
pub mod reach;
pub mod weight;

pub const KNN_ASSEMBLY: &str = "knn_assembly";
pub const KNN_POOL: &str = "knn_pool";
pub const KNN_SPLIT: &str = "knn_split";
pub const KNN_ATTACH: &str = "knn_attach";

pub type Graph = sprs::CsMatI<f32, u32, usize>;

pub(crate) fn row_of(graph: &Graph, row: usize) -> (&[u32], &[f32]) {
    let start = graph.indptr().index(row);
    let end = graph.indptr().index(row + 1);
    (&graph.indices()[start..end], &graph.data()[start..end])
}

pub(crate) fn transpose<T: Copy + Default>(
    graph: &Graph,
    value: impl Fn(usize, usize) -> T,
) -> (Vec<usize>, Vec<T>) {
    let nodes = graph.rows();
    let mut starts = vec![0usize; nodes + 1];
    for column in graph.indices() {
        starts[*column as usize + 1] += 1;
    }
    for at in 0..nodes {
        starts[at + 1] += starts[at];
    }
    let mut next = starts.clone();
    let mut values = vec![T::default(); graph.indices().len()];
    for row in 0..nodes {
        let start = graph.indptr().index(row);
        for (offset, column) in row_of(graph, row).0.iter().enumerate() {
            values[next[*column as usize]] = value(row, start + offset);
            next[*column as usize] += 1;
        }
    }
    (starts, values)
}

// `Level::from_graph` and `label_propagation` size their vectors by `graph.rows()` and index
// them by node id, so a slice has to come back renumbered rather than masked.
pub fn induced(graph: &Graph, indices: &[usize]) -> Graph {
    let mut indptr = Vec::with_capacity(indices.len() + 1);
    let (mut columns, mut data) = (Vec::new(), Vec::new());
    indptr.push(0);
    for node in indices {
        let (held, weights) = row_of(graph, *node);
        for (column, weight) in held.iter().zip(weights) {
            if let Ok(position) = indices.binary_search(&(*column as usize)) {
                columns.push(position as u32);
                data.push(*weight);
            }
        }
        indptr.push(columns.len());
    }
    Graph::new((indices.len(), indices.len()), indptr, columns, data)
}

// The assembler's own adjacency is weak on its own, a coin flip as a must-link on the one real
// set with a gold, so it joins the neighbour graph as another edge rather than as a constraint.
pub fn linked(
    graph: Graph,
    links: &[crate::assembly_graph::Link],
    indices: &[usize],
    weight: f32,
) -> Graph {
    let mut wanted = links
        .iter()
        .filter_map(|link| {
            let from = indices.binary_search(&link.from).ok()? as u32;
            let to = indices.binary_search(&link.to).ok()? as u32;
            (from != to).then_some([(from, to), (to, from)])
        })
        .flatten()
        .collect::<Vec<_>>();
    if wanted.is_empty() {
        return graph;
    }
    wanted.sort_unstable();
    wanted.dedup();
    // The neighbour graph is symmetric, so a pair missing from one row is missing from both.
    let rows = graph.rows();
    let mut indptr = Vec::with_capacity(rows + 1);
    let mut out_columns = Vec::with_capacity(graph.nnz() + wanted.len());
    let mut out_data = Vec::with_capacity(graph.nnz() + wanted.len());
    indptr.push(0);
    let mut next = wanted.into_iter().peekable();
    for row in 0..rows {
        let (columns, weights) = row_of(&graph, row);
        let row = row as u32;
        for (column, edge) in columns.iter().zip(weights) {
            while let Some((_, added)) = next.next_if(|pair| *pair < (row, *column)) {
                out_columns.push(added);
                out_data.push(weight);
            }
            let raised = next.next_if_eq(&(row, *column)).is_some();
            out_columns.push(*column);
            out_data.push(if raised { edge.max(weight) } else { *edge });
        }
        while let Some((_, added)) = next.next_if(|pair| pair.0 == row) {
            out_columns.push(added);
            out_data.push(weight);
        }
        indptr.push(out_columns.len());
    }
    Graph::new((rows, rows), indptr, out_columns, out_data)
}
