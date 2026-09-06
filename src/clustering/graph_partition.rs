use std::borrow::Cow;

use log::info;
use rand::{SeedableRng, rngs::StdRng, seq::SliceRandom};

use crate::embedding::{Graph, intersect::row_of};

const MAX_ROUNDS: usize = 50;

/// HDBSCAN reads the layout, the other two read the graph it was built from. Neither graph
/// source has a noise label, so the eject is what refuses a contig under them.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Partition {
    Hdbscan,
    LabelProp,
    Leiden,
    Infomap,
    Sbm,
    #[default]
    Auto,
}

pub const PARTITION_NAMES: [&str; 6] = ["auto", "hdbscan", "labelprop", "leiden", "infomap", "sbm"];

/// No contig length statistic separates the two arms: CAMI I medium and low agree on median
/// and want opposite ones. Total assembly does, over the eight datasets measured.
const LARGE_ASSEMBLY_BP: usize = 350_000_000;

impl Partition {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "auto" => Some(Self::Auto),
            "hdbscan" => Some(Self::Hdbscan),
            "labelprop" => Some(Self::LabelProp),
            "leiden" => Some(Self::Leiden),
            "infomap" => Some(Self::Infomap),
            "sbm" => Some(Self::Sbm),
            _ => None,
        }
    }

    pub fn resolve(self, lengths: &[usize]) -> Self {
        if self != Self::Auto {
            return self;
        }
        let total = lengths.iter().sum::<usize>();
        let chosen = if total >= LARGE_ASSEMBLY_BP {
            Self::LabelProp
        } else {
            Self::Leiden
        };
        info!(
            "Assembly {} Mbp past the filter, partitioning with {chosen:?}",
            total / 1_000_000
        );
        chosen
    }

    pub fn reads_graph(&self) -> bool {
        *self != Self::Hdbscan
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum NodeSize {
    Count,
    #[default]
    Bp,
}

pub const NODE_SIZE_NAMES: [&str; 2] = ["count", "bp"];

pub struct SizedGraph<'a> {
    pub graph: Cow<'a, Graph>,
    pub sizes: Option<Vec<f64>>,
}

impl NodeSize {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "count" => Some(Self::Count),
            "bp" => Some(Self::Bp),
            _ => None,
        }
    }

    pub fn apply<'a>(self, graph: &'a Graph, lengths: &[usize]) -> SizedGraph<'a> {
        match self {
            Self::Count => SizedGraph {
                graph: Cow::Borrowed(graph),
                sizes: None,
            },
            Self::Bp => SizedGraph {
                graph: Cow::Owned(bp_weighted(graph, lengths)),
                sizes: Some(lengths.iter().map(|length| *length as f64).collect()),
            },
        }
    }
}

/// Mass in bases alone shreds a genome held in a few large contigs, and edges scaled by each
/// end's own length hand a closed contig its whole length of foreign pull. The geometric mean
/// is the one scaling that keeps a genome-sized contig out of its neighbours' bin without either.
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

pub struct Weights {
    pub degrees: Vec<f64>,
    pub total: f64,
}

impl Weights {
    pub fn of(graph: &Graph) -> Self {
        let degrees = (0..graph.rows())
            .map(|row| row_of(graph, row).1.iter().map(|w| *w as f64).sum())
            .collect::<Vec<f64>>();
        let total = degrees.iter().sum::<f64>() / 2.0;
        Self { degrees, total }
    }
}

/// Self loops are dropped so the degree convention matches `label_propagation` and
/// `Level::from_graph`, which both skip them.
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

pub fn visit_order(nodes: usize, seed: u64) -> Vec<usize> {
    let mut order = (0..nodes).collect::<Vec<_>>();
    order.shuffle(&mut StdRng::seed_from_u64(seed));
    order
}

pub fn label_propagation(graph: &Graph, seed: u64) -> Vec<i32> {
    let nodes = graph.rows();
    let mut labels = (0..nodes as i32).collect::<Vec<i32>>();
    let order = visit_order(nodes, seed);
    let mut incident = Incident::new(nodes);

    for _ in 0..MAX_ROUNDS {
        let mut moved = false;
        for node in &order {
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
            }
        }
        if !moved {
            break;
        }
    }

    compact(&labels)
}

/// Reused across visits: gathering runs once per node per round here and once per queue pop
/// in Leiden, so the map it replaces was built tens of millions of times a run.
pub struct Incident {
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

pub fn compact(labels: &[i32]) -> Vec<i32> {
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
