use std::collections::{HashMap, HashSet};

use crate::clustering::clusterer::HDBSCANResult;
use crate::embedding::{Graph, intersect};

/// Positions a must-link chain joins partition as one node, so their mass counts once and whole
/// rather than twice at half, which is the only way the constraint can change a bp weighted move.
pub struct Contraction {
    of: Vec<usize>,
    members: Vec<Vec<usize>>,
}

impl Contraction {
    pub fn new(component: &[usize], contigs: &[usize]) -> Option<Self> {
        let mut group_of: HashMap<usize, usize> = HashMap::new();
        let mut members: Vec<Vec<usize>> = Vec::new();
        let mut of = Vec::with_capacity(contigs.len());
        for (position, contig) in contigs.iter().enumerate() {
            let group = *group_of.entry(component[*contig]).or_insert_with(|| {
                members.push(Vec::new());
                members.len() - 1
            });
            members[group].push(position);
            of.push(group);
        }
        (members.len() < contigs.len()).then_some(Self { of, members })
    }

    pub fn len(&self) -> usize {
        self.members.len()
    }

    pub fn graph(&self, graph: &Graph) -> Graph {
        let groups = self.members.len();
        let mut indptr = Vec::with_capacity(groups + 1);
        let mut indices = Vec::with_capacity(graph.nnz());
        let mut data = Vec::with_capacity(graph.nnz());
        indptr.push(0);
        let mut weights: HashMap<u32, f32> = HashMap::new();
        for (group, positions) in self.members.iter().enumerate() {
            weights.clear();
            for position in positions {
                let (targets, found) = intersect::row_of(graph, *position);
                for (target, weight) in targets.iter().zip(found) {
                    let other = self.of[*target as usize];
                    if other != group {
                        *weights.entry(other as u32).or_default() += *weight;
                    }
                }
            }
            let mut row = weights
                .iter()
                .map(|(at, weight)| (*at, *weight))
                .collect::<Vec<_>>();
            row.sort_unstable_by_key(|(at, _)| *at);
            indices.extend(row.iter().map(|(at, _)| *at));
            data.extend(row.iter().map(|(_, weight)| *weight));
            indptr.push(indices.len());
        }
        sprs::CsMatI::new((groups, groups), indptr, indices, data)
    }

    pub fn lengths(&self, lengths: &[usize]) -> Vec<usize> {
        self.members
            .iter()
            .map(|positions| positions.iter().map(|at| lengths[*at]).sum())
            .collect()
    }

    pub fn expand(&self, result: HDBSCANResult) -> HDBSCANResult {
        let spread = |groups: &HashSet<usize>| {
            groups
                .iter()
                .flat_map(|group| self.members[*group].iter().copied())
                .collect::<HashSet<_>>()
        };
        HDBSCANResult {
            cluster_map: result
                .cluster_map
                .iter()
                .map(|(label, groups)| (*label, spread(groups)))
                .collect(),
            outliers: spread(&result.outliers),
            score: result.score,
        }
    }
}
