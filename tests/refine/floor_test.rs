//! Whether a bin holding more than one genome-sized contig comes apart into them, and
//! where genome-sized is read from.

use std::collections::BTreeMap;

use ndarray::Array2;
use rosella::clustering::graph_partition::{NodeSize, Partition};
use rosella::clustering::objective::ObjectiveChoice;
use rosella::embedding::features::ContigFeatures;
use rosella::embedding::umap::EmbedOverrides;
use rosella::refine::bin_stats::LevelSource;
use rosella::refine::floor::floor;
use rosella::refine::splitter::{RefineSettings, Refiner};
use rosella::seeds::Seeds;

const MIN_BIN_SIZE: usize = 200_000;

fn features(lengths: &[usize]) -> (Array2<f64>, Array2<f64>) {
    (
        Array2::zeros((lengths.len(), 2)),
        Array2::zeros((lengths.len(), 6)),
    )
}

/// A lone short contig is an outlier, not a genome, so it says nothing about genome size.
#[test]
fn the_floor_is_half_the_median_closed_genome() {
    let lengths = [3_000_000, 2_000_000, 50_000, 2_600_000, 20_000, 10_000];
    let (coverage, tnf) = features(&lengths);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let bins = BTreeMap::from([(0, vec![0]), (1, vec![1]), (2, vec![2]), (3, vec![4, 5])]);

    assert_eq!(floor(&features, &bins, &[3], MIN_BIN_SIZE), Some(1_300_000));
    assert_eq!(
        floor(&features, &BTreeMap::new(), &[2], MIN_BIN_SIZE,),
        None
    );
}

/// A closed genome that collected short contigs still carries most of its bin, and a bin of

/// recover turns refinement off by passing 0 rounds, and the pool reads the floor the refiner
/// measured, so returning early without measuring it made `--no-refine` two ablations at once.
#[test]
fn the_floor_is_measured_even_when_no_round_runs() {
    let lengths = [3_000_000, 2_000_000, 2_600_000, 10_000];
    let (coverage, tnf) = features(&lengths);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let objective = ObjectiveChoice::parse("codelength").unwrap().build();
    let settings = RefineSettings {
        min_bin_size: MIN_BIN_SIZE,
        max_bin_size: 15_000_000,
        n_neighbours: 100,
        max_retries: 0,
        seeds: Seeds {
            knn: 42,
            sample: 42,
            partition: 42,
        },
        max_contamination: None,
        overrides: EmbedOverrides::default(),
        bisect: false,
        levels: LevelSource::Derived,
        level_quantile: 0.75,
        partition: Partition::Auto.resolve(&lengths),
        node_size: NodeSize::Count,
        partition_resolution: None,
        partition_theta: None,
    };
    let bins = BTreeMap::from([(0usize, vec![0]), (1, vec![1]), (2, vec![2]), (3, vec![3])]);
    let mut refiner = Refiner::new(features, &objective, settings, bins, Vec::new());

    assert_eq!(refiner.run(), 0);
    assert_eq!(refiner.genome_floor, Some(1_300_000));
}
