use std::collections::{HashMap, HashSet};

use ndarray::array;
use rosella::clustering::clusterer::Partitioning;
use rosella::clustering::graph_partition::Partition;
use rosella::coverage::coverage_table::CoverageTable;
use rosella::kmers::kmer_counting::{KmerFrequencyTable, KmerSizes};
use rosella::recover::recover_engine::attract::Attractors;

fn tables(lengths: &[usize]) -> (CoverageTable, KmerFrequencyTable) {
    let names = (0..lengths.len()).map(|at| format!("c{at}")).collect::<Vec<_>>();
    let rows = lengths.len();
    let coverage = CoverageTable {
        table: ndarray::Array2::from_shape_fn((rows, 2), |(row, _)| row as f64),
        average_depths: vec![1.0; rows],
        contig_names: names.clone(),
        contig_lengths: lengths.to_vec(),
        sample_names: vec!["s".to_string()],
    };
    let composition = KmerFrequencyTable::new(
        "4".parse::<KmerSizes>().unwrap(),
        ndarray::Array2::from_shape_fn((rows, 1), |(row, _)| row as f64),
        names,
    );
    (coverage, composition)
}

#[test]
fn short_rows_leave_the_tables_and_every_bin() {
    let (mut coverage, mut composition) = tables(&[2000, 800, 1600, 900]);
    let view = Attractors::split(&mut coverage, &mut composition, 1500, false)
        .unwrap()
        .unwrap();
    assert_eq!(coverage.contig_lengths, vec![2000, 1600]);
    assert_eq!(composition.kmer_table, array![[0.0], [2.0]]);

    let held = view.binnable(Partitioning {
        cluster_map: HashMap::from([
            (0, HashSet::from([0, 1])),
            (1, HashSet::from([3])),
            (2, HashSet::from([2])),
        ]),
        outliers: HashSet::new(),
        score: None,
        arm: Partition::Leiden,
        seed: 0,
    });
    let mut bins = held
        .cluster_map
        .into_values()
        .map(|members| {
            let mut members = members.into_iter().collect::<Vec<_>>();
            members.sort_unstable();
            members
        })
        .collect::<Vec<_>>();
    bins.sort_unstable();
    assert_eq!(bins, vec![vec![0], vec![1]]);
}

#[test]
fn nothing_under_the_cutoff_leaves_no_view() {
    let (mut coverage, mut composition) = tables(&[2000, 1600]);
    assert!(
        Attractors::split(&mut coverage, &mut composition, 1500, false)
            .unwrap()
            .is_none()
    );
}
