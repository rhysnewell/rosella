use crate::embedding::knn::KnnGraph;

/// Two candidates a hundredth of a bin apart are the same proposal to the bar, so a lineage
/// only re-enters the heap once it has grown enough to be a different answer.
const GROWTH: f64 = 1.01;

struct Forest {
    parent: Vec<usize>,
    members: Vec<Vec<usize>>,
    bases: Vec<usize>,
    emitted: Vec<usize>,
}

impl Forest {
    fn new(order: &[usize], length: &impl Fn(usize) -> usize) -> Self {
        Self {
            parent: (0..order.len()).collect(),
            members: order.iter().map(|contig| vec![*contig]).collect(),
            bases: order.iter().map(|contig| length(*contig)).collect(),
            emitted: vec![0; order.len()],
        }
    }

    fn root(&mut self, node: usize) -> usize {
        let mut node = node;
        while self.parent[node] != node {
            self.parent[node] = self.parent[self.parent[node]];
            node = self.parent[node];
        }
        node
    }

    fn union(&mut self, left: usize, right: usize) -> Option<usize> {
        let (left, right) = (self.root(left), self.root(right));
        if left == right {
            return None;
        }
        let (keep, drop) = match self.members[left].len() >= self.members[right].len() {
            true => (left, right),
            false => (right, left),
        };
        let taken = std::mem::take(&mut self.members[drop]);
        self.members[keep].extend(taken);
        self.bases[keep] += self.bases[drop];
        self.emitted[keep] = self.emitted[keep].max(self.emitted[drop]);
        self.parent[drop] = keep;
        Some(keep)
    }
}

fn edges(knn: &KnnGraph) -> Vec<(f32, u32, u32)> {
    let mut held = Vec::with_capacity(knn.indices.len());
    for (row, neighbours) in knn.indices.rows().into_iter().enumerate() {
        for (column, neighbour) in neighbours.iter().enumerate() {
            let (left, right) = match (row as u32, *neighbour) {
                (row, neighbour) if row == neighbour => continue,
                (row, neighbour) if row < neighbour => (row, neighbour),
                (row, neighbour) => (neighbour, row),
            };
            held.push((knn.dists[[row, column]], left, right));
        }
    }
    held.sort_unstable_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then_with(|| (left.1, left.2).cmp(&(right.1, right.2)))
    });
    held
}

/// A genome holding a handful of nodes in a dense graph is never a community at any rung, but
/// it is a whole subtree of the same graph's merge order.
pub fn candidates(
    knn: &KnnGraph,
    order: &[usize],
    length: impl Fn(usize) -> usize,
    floor: usize,
    ceiling: usize,
) -> Vec<Vec<usize>> {
    if order.len() != knn.n_points() {
        return Vec::new();
    }
    let mut forest = Forest::new(order, &length);
    let mut found = Vec::new();
    for (_, left, right) in edges(knn) {
        let Some(root) = forest.union(left as usize, right as usize) else {
            continue;
        };
        let bases = forest.bases[root];
        if bases < floor || bases > ceiling {
            continue;
        }
        if (bases as f64) < forest.emitted[root] as f64 * GROWTH {
            continue;
        }
        forest.emitted[root] = bases;
        let mut group = forest.members[root].clone();
        group.sort_unstable();
        found.push(group);
    }
    found
}
