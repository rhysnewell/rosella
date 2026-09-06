//! Label propagation and Leiden both read the fuzzy simplicial set the layout was built from,
//! so what they have to recover is the block structure a k-NN graph over clustered contigs
//! carries. The fixture plants that structure so a partition can be checked against a known
//! answer rather than against another partition.

use rosella::clustering::graph_partition::label_propagation;
use rosella::clustering::infomap::infomap;
use rosella::clustering::leiden::{leiden, resolutions};
use rosella::clustering::sbm::sbm;
use sprs::{CsMatI, TriMatI};
use std::collections::{HashMap, HashSet};

const PER_BLOCK: usize = 40;
const BLOCKS: usize = 12;

fn blocked_graph(blocks: usize, bridge: f32) -> CsMatI<f32, u32, usize> {
    let n = blocks * PER_BLOCK;
    let mut triplets = TriMatI::new((n, n));
    let mut add = |i: usize, j: usize, weight: f32| {
        triplets.add_triplet(i, j, weight);
        triplets.add_triplet(j, i, weight);
    };
    for block in 0..blocks {
        let start = block * PER_BLOCK;
        for i in start..start + PER_BLOCK {
            for j in i + 1..start + PER_BLOCK {
                add(i, j, 1.0);
            }
        }
        add(start, (start + PER_BLOCK) % n, bridge);
    }
    triplets.to_csr()
}

fn planted(labels: &[i32]) -> f64 {
    let mut agree = 0.0;
    let mut total = 0.0;
    for i in 0..labels.len() {
        for j in i + 1..labels.len() {
            let same_block = i / PER_BLOCK == j / PER_BLOCK;
            if same_block == (labels[i] == labels[j]) {
                agree += 1.0;
            }
            total += 1.0;
        }
    }
    agree / total
}

fn communities(labels: &[i32]) -> usize {
    labels.iter().collect::<HashSet<_>>().len()
}

#[test]
fn label_propagation_recovers_the_planted_blocks() {
    let graph = blocked_graph(BLOCKS, 0.01);
    let labels = label_propagation(&graph, 42);
    let agreement = planted(&labels);
    assert!(
        agreement > 0.99,
        "label propagation agreed with the planted blocks on only {agreement:.4} of pairs, \
         over {} communities against {BLOCKS} planted",
        communities(&labels)
    );
}

#[test]
fn leiden_recovers_the_planted_blocks() {
    let graph = blocked_graph(BLOCKS, 0.01);
    let gamma = resolutions(&graph, None, 8);
    let best = gamma
        .iter()
        .map(|resolution| {
            let labels = leiden(&graph, None, *resolution, None, 42);
            (planted(&labels), communities(&labels))
        })
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .expect("the resolution ladder is never empty");
    assert!(
        best.0 > 0.99,
        "the best of {} resolutions agreed on only {:.4} of pairs, over {} communities \
         against {BLOCKS} planted",
        gamma.len(),
        best.0,
        best.1
    );
}

#[test]
fn a_single_blob_yields_one_community() {
    let graph = blocked_graph(1, 0.0);
    assert_eq!(communities(&label_propagation(&graph, 42)), 1);
}

#[test]
fn every_node_is_assigned() {
    let graph = blocked_graph(BLOCKS, 0.01);
    let labels = leiden(&graph, None, resolutions(&graph, None, 8)[4], None, 42);
    assert_eq!(labels.len(), BLOCKS * PER_BLOCK);
    assert!(
        labels.iter().all(|label| *label >= 0),
        "a graph partition has no noise label, so nothing may come back negative"
    );
}

#[test]
fn resolution_trades_community_count_against_size() {
    let graph = blocked_graph(BLOCKS, 0.01);
    let ladder = resolutions(&graph, None, 6);
    let counts = ladder
        .iter()
        .map(|resolution| communities(&leiden(&graph, None, *resolution, None, 42)))
        .collect::<Vec<_>>();
    assert!(
        counts.first() <= counts.last(),
        "raising the resolution should not cut the community count: {counts:?} over {ladder:?}"
    );
}

#[test]
fn a_bridged_pair_stays_apart_at_the_resolution_that_splits_it() {
    let graph = blocked_graph(2, 0.5);
    let ladder = resolutions(&graph, None, 8);
    let splits = ladder
        .iter()
        .any(|resolution| communities(&leiden(&graph, None, *resolution, None, 42)) == 2);
    assert!(
        splits,
        "no resolution on the ladder parted two cliques joined by a single weak edge"
    );
}

#[test]
fn label_propagation_is_reproducible_across_a_seed_change() {
    let graph = blocked_graph(BLOCKS, 0.01);
    let a = label_propagation(&graph, 42);
    let b = label_propagation(&graph, 7);
    let pairs = |labels: &[i32]| {
        let mut groups: HashMap<i32, Vec<usize>> = HashMap::new();
        for (node, label) in labels.iter().enumerate() {
            groups.entry(*label).or_default().push(node);
        }
        let mut sets = groups.into_values().collect::<Vec<_>>();
        sets.sort();
        sets
    };
    assert_eq!(
        pairs(&a),
        pairs(&b),
        "well separated blocks should not depend on the visit order"
    );
}

/// Infomap optimises the map equation directly rather than ranking a CPM ladder on it, so the
/// question is whether the search that needs no resolution still finds the planted answer.
#[test]
fn infomap_recovers_the_planted_blocks() {
    let graph = blocked_graph(BLOCKS, 0.01);
    let labels = infomap(&graph, 42);
    let agreement = planted(&labels);
    assert!(
        agreement > 0.99 && communities(&labels) == BLOCKS,
        "infomap agreed with the planted blocks on only {agreement:.4} of pairs, \
         over {} communities against {BLOCKS} planted",
        communities(&labels)
    );
}

/// A graph with no block structure has nothing to describe more cheaply than one module, and
/// the map equation is what says so. Leiden needs the ladder selector to reach the same answer.
#[test]
fn infomap_leaves_a_single_blob_whole() {
    let graph = blocked_graph(1, 0.0);
    assert_eq!(communities(&infomap(&graph, 42)), 1);
}

/// The block model selects its own block count under minimum description length, so unlike the
/// ladder arms it has to land on the planted number without being told what to look for.
#[test]
fn the_block_model_recovers_the_planted_blocks() {
    let graph = blocked_graph(BLOCKS, 0.01);
    let labels = sbm(&graph, 42);
    let agreement = planted(&labels);
    assert!(
        agreement > 0.99 && communities(&labels) == BLOCKS,
        "the block model agreed with the planted blocks on only {agreement:.4} of pairs, \
         over {} communities against {BLOCKS} planted",
        communities(&labels)
    );
}
