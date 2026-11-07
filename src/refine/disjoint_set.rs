pub struct DisjointSet {
    parent: Vec<usize>,
}

impl DisjointSet {
    pub fn new(nodes: usize) -> Self {
        Self {
            parent: (0..nodes).collect(),
        }
    }

    pub fn find(&mut self, mut node: usize) -> usize {
        while self.parent[node] != node {
            self.parent[node] = self.parent[self.parent[node]];
            node = self.parent[node];
        }
        node
    }

    pub fn attach(&mut self, root: usize, under: usize) {
        self.parent[root] = under;
    }
}
