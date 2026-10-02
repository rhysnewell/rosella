use rand::{SeedableRng, rngs::StdRng, seq::SliceRandom};

use crate::embedding::{Graph, row_of};

const MAX_ROUNDS: usize = 50;

// Neither source has a noise label, so the eject is what refuses a contig under them.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum Partition {
    #[value(name = "labelprop")]
    LabelProp,
    Leiden,
    #[default]
    Both,
}

impl Partition {
    // Label propagation returns one labelling, so it is also the arm that offers no ladder.
    pub fn runs_leiden(&self) -> bool {
        *self != Self::LabelProp
    }

    pub fn runs_labelprop(&self) -> bool {
        matches!(self, Self::LabelProp | Self::Both)
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::LabelProp => "labelprop",
            Self::Leiden => "leiden",
            Self::Both => "both",
        }
    }

    // The splitter cuts one bin at a time and keeps whatever codelength ranks first. Running
    // the second arm here was measured over 13 sets and moved nothing, cami_i_high included.
    pub fn for_split(self) -> Self {
        match self {
            Self::Both => Self::Leiden,
            chosen => chosen,
        }
    }
}

pub(super) struct SizedGraph {
    pub graph: Graph,
    pub sizes: Vec<f64>,
}

pub(super) fn sized(graph: &Graph, lengths: &[usize]) -> SizedGraph {
    SizedGraph {
        graph: bp_weighted(graph, lengths),
        sizes: lengths.iter().map(|length| *length as f64).collect(),
    }
}

// Bases alone shred a genome held in a few large contigs, and one end's own length gives a closed
// contig its whole length of foreign pull. The geometric mean does neither.
fn bp_weighted(graph: &Graph, lengths: &[usize]) -> Graph {
    let boundaries = (0..graph.rows())
        .map(|row| (graph.indptr().index(row), graph.indptr().index(row + 1)))
        .collect::<Vec<_>>();
    let targets = graph.indices();
    let mut weighted = graph.clone();
    let data = weighted.data_mut();
    for (row, (start, end)) in boundaries.into_iter().enumerate() {
        let own = lengths[row] as f32;
        for (slot, target) in data[start..end].iter_mut().zip(&targets[start..end]) {
            *slot *= (own * lengths[*target as usize] as f32).sqrt();
        }
    }
    weighted
}

// Self loops are kept here, unlike `node_degrees`, because the ladder spans a range of
// community sizes rather than scoring a partition.
pub(super) fn edge_weight_total(graph: &Graph) -> f64 {
    (0..graph.rows())
        .map(|row| row_of(graph, row).1.iter().map(|w| *w as f64).sum::<f64>())
        .sum::<f64>()
        / 2.0
}

// Self loops are dropped so the degree convention matches `label_propagation` and
// `Level::from_graph`, which both skip them.
pub fn node_degrees(graph: &Graph) -> Vec<f64> {
    (0..graph.rows())
        .map(|row| {
            let (targets, weights) = row_of(graph, row);
            targets
                .iter()
                .zip(weights)
                .filter(|(target, _)| **target as usize != row)
                .map(|(_, weight)| *weight as f64)
                .sum()
        })
        .collect()
}

pub(super) fn visit_order(nodes: usize, seed: u64) -> Vec<usize> {
    let mut order = (0..nodes).collect::<Vec<_>>();
    order.shuffle(&mut StdRng::seed_from_u64(seed));
    order
}

// A node whose neighbours kept their labels since its last visit would pick the label it holds,
// so only a node one of them left is visited. The labelling is the one visiting every node gives.
pub fn label_propagation(graph: &Graph, seed: u64) -> Vec<i32> {
    let nodes = graph.rows();
    let mut labels = (0..nodes as i32).collect::<Vec<i32>>();
    let order = visit_order(nodes, seed);
    let mut incident = Incident::new(nodes);
    let (starts, readers) = readers(graph);
    let mut stale = vec![true; nodes];

    for _ in 0..MAX_ROUNDS {
        let mut moved = false;
        for node in &order {
            if !std::mem::take(&mut stale[*node]) {
                continue;
            }
            let (neighbours, weights) = row_of(graph, *node);
            if neighbours.is_empty() {
                continue;
            }
            incident.gather(
                neighbours
                    .iter()
                    .zip(weights)
                    .filter(|(target, _)| **target as usize != *node)
                    .map(|(target, weight)| (labels[*target as usize] as usize, *weight as f64)),
            );
            let best = incident
                .iter()
                .max_by(|a, b| a.1.total_cmp(&b.1).then_with(|| b.0.cmp(&a.0)));
            if let Some((label, _)) = best
                && labels[*node] != label as i32
            {
                labels[*node] = label as i32;
                moved = true;
                for reader in &readers[starts[*node]..starts[*node + 1]] {
                    stale[*reader as usize] = true;
                }
            }
        }
        if !moved {
            break;
        }
    }

    compact(&labels)
}

// Taken from the rows rather than assumed symmetric, so a move marks every node that reads it.
fn readers(graph: &Graph) -> (Vec<usize>, Vec<u32>) {
    let nodes = graph.rows();
    let mut starts = vec![0usize; nodes + 1];
    for column in graph.indices() {
        starts[*column as usize + 1] += 1;
    }
    for at in 0..nodes {
        starts[at + 1] += starts[at];
    }
    let mut next = starts.clone();
    let mut readers = vec![0u32; graph.indices().len()];
    for row in 0..nodes {
        for column in row_of(graph, row).0 {
            readers[next[*column as usize]] = row as u32;
            next[*column as usize] += 1;
        }
    }
    (starts, readers)
}

// Reused across visits: gathering runs once per node per round here and once per queue pop
// in Leiden, so the map it replaces was built tens of millions of times a run.
pub(super) struct Incident {
    totals: Vec<f64>,
    stamp: Vec<u64>,
    touched: Vec<usize>,
    generation: u64,
}

impl Incident {
    pub fn new(communities: usize) -> Self {
        Self {
            totals: vec![0.0; communities],
            stamp: vec![0; communities],
            touched: Vec::new(),
            generation: 0,
        }
    }

    pub fn gather(&mut self, incident: impl Iterator<Item = (usize, f64)>) {
        self.generation += 1;
        self.touched.clear();
        for (community, weight) in incident {
            if self.stamp[community] != self.generation {
                self.stamp[community] = self.generation;
                self.totals[community] = 0.0;
                self.touched.push(community);
            }
            self.totals[community] += weight;
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, f64)> + '_ {
        self.touched
            .iter()
            .map(move |community| (*community, self.totals[*community]))
    }

    pub fn get(&self, community: usize) -> f64 {
        if self.stamp[community] == self.generation {
            self.totals[community]
        } else {
            0.0
        }
    }
}

pub(super) fn compact(labels: &[i32]) -> Vec<i32> {
    let Some(highest) = labels.iter().copied().max() else {
        return Vec::new();
    };
    let lowest = labels.iter().copied().min().unwrap_or(highest);
    let mut first = vec![-1i32; (highest - lowest) as usize + 1];
    let mut next = 0i32;
    labels
        .iter()
        .map(|label| {
            let slot = &mut first[(*label - lowest) as usize];
            if *slot < 0 {
                *slot = next;
                next += 1;
            }
            *slot
        })
        .collect()
}
