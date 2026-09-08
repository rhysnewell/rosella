use std::collections::{HashMap, HashSet};

use rosella::clustering::clusterer::HDBSCANResult;
use rosella::clustering::contract::Contraction;
use rosella::embedding::Graph;
use sprs::TriMatI;

fn graph(edges: &[(usize, usize, f32)], nodes: usize) -> Graph {
    let mut triplets = TriMatI::new((nodes, nodes));
    for (one, other, weight) in edges {
        triplets.add_triplet(*one, *other, *weight);
        triplets.add_triplet(*other, *one, *weight);
    }
    triplets.to_csr()
}

fn row(graph: &Graph, at: usize) -> Vec<(usize, f32)> {
    let view = graph.outer_view(at).unwrap();
    view.iter()
        .map(|(index, weight)| (index, *weight))
        .collect()
}

#[test]
fn a_component_of_one_contig_each_contracts_nothing() {
    assert!(Contraction::new(&[0, 1, 2], &[0, 1, 2]).is_none());
}

#[test]
fn joined_positions_become_one_node_that_carries_both_edges() {
    let contraction = Contraction::new(&[0, 0, 2], &[0, 1, 2]).unwrap();
    let held = contraction.graph(&graph(&[(0, 1, 5.0), (0, 2, 1.0), (1, 2, 2.0)], 3));

    assert_eq!(contraction.len(), 2);
    assert_eq!(held.rows(), 2);
    assert_eq!(row(&held, 0), vec![(1usize, 3.0)]);
    assert_eq!(contraction.lengths(&[100, 250, 400]), vec![350, 400]);
}

#[test]
fn expanding_puts_every_joined_position_in_the_same_cluster() {
    let contraction = Contraction::new(&[0, 0, 2], &[0, 1, 2]).unwrap();
    let result = HDBSCANResult {
        cluster_map: HashMap::from([(0, HashSet::from([0]))]),
        outliers: HashSet::from([1]),
        score: 0.5,
    };

    let expanded = contraction.expand(result);

    assert_eq!(expanded.cluster_map[&0], HashSet::from([0, 1]));
    assert_eq!(expanded.outliers, HashSet::from([2]));
}

#[test]
fn contigs_outside_the_subset_do_not_shift_the_grouping() {
    let contraction = Contraction::new(&[0, 0, 2, 2], &[1, 2, 3]).unwrap();

    assert_eq!(contraction.len(), 2);
    assert_eq!(contraction.lengths(&[10, 20, 30]), vec![10, 50]);
}
