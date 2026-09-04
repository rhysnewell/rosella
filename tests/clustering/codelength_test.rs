//! Modularity picks the coarsest rung of the ladder because of its resolution limit, so the
//! fixture is the ring of cliques that limit is defined on. With 24 cliques of 4 the merged
//! partition beats the planted one under modularity, and the point of the map equation is
//! that it does not.

use rosella::clustering::codelength::codelength_saving;
use rosella::clustering::modularity::modularity;
use sprs::{CsMatI, TriMatI};

const CLIQUES: usize = 24;
const SIZE: usize = 4;
const NODES: usize = CLIQUES * SIZE;

fn ring_of_cliques() -> CsMatI<f32, u32, usize> {
    let mut triplets = TriMatI::new((NODES, NODES));
    let mut add = |i: usize, j: usize| {
        triplets.add_triplet(i, j, 1.0);
        triplets.add_triplet(j, i, 1.0);
    };
    for clique in 0..CLIQUES {
        let start = clique * SIZE;
        for i in start..start + SIZE {
            for j in i + 1..start + SIZE {
                add(i, j);
            }
        }
        add(start, (start + SIZE) % NODES);
    }
    triplets.to_csr()
}

fn grouped(per_label: usize) -> Vec<i32> {
    (0..NODES).map(|node| (node / per_label) as i32).collect()
}

#[test]
fn a_single_community_saves_nothing() {
    let saving = codelength_saving(&ring_of_cliques(), &vec![0i32; NODES]);
    assert!(
        saving.abs() < 1e-12,
        "one community should save exactly nothing, got {saving}"
    );
}

#[test]
fn the_saving_holds_the_cliques_modularity_merges() {
    let graph = ring_of_cliques();
    let planted = grouped(SIZE);
    let merged = grouped(SIZE * 2);

    let merged_modularity = modularity(&graph, &merged, 1.0);
    let planted_modularity = modularity(&graph, &planted, 1.0);
    assert!(
        merged_modularity > planted_modularity,
        "the fixture is meant to trip modularity's resolution limit, but planted \
         {planted_modularity} already beats merged {merged_modularity}"
    );

    let planted_saving = codelength_saving(&graph, &planted);
    let merged_saving = codelength_saving(&graph, &merged);
    assert!(
        planted_saving > merged_saving,
        "planted {planted_saving} should beat merged {merged_saving}"
    );
}

#[test]
fn the_saving_refuses_to_split_a_clique() {
    let graph = ring_of_cliques();
    let planted = codelength_saving(&graph, &grouped(SIZE));
    let split = codelength_saving(&graph, &grouped(SIZE / 2));
    assert!(
        planted > split,
        "planted {planted} should beat split {split}"
    );
}
