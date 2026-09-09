//! Whether the re-clustering that came back is worth taking. Whether the bin was worth
//! re-clustering in the first place is `bar_test`.

use std::collections::BTreeMap;

use ndarray::Array2;
use rosella::clustering::graph_partition::{NodeSize, Partition};
use rosella::clustering::objective::ObjectiveChoice;
use rosella::embedding::features::ContigFeatures;
use rosella::embedding::umap::EmbedOverrides;
use rosella::refine::bin_stats::LevelSource;
use rosella::refine::gates::SplitRejection;
use rosella::refine::proposal::{judge_split, leaves_two_standing};
use rosella::refine::splitter::{RefineSettings, Refiner};
use rosella::seeds::Seeds;

const CONTIG_LENGTH: usize = 100_000;

fn size_of(cluster: &[usize]) -> usize {
    cluster.len() * CONTIG_LENGTH
}

fn cluster(size: usize, offset: usize) -> Vec<usize> {
    (offset..offset + size).collect()
}

#[test]
fn split_rejections() {
    let cases: [(Vec<Vec<usize>>, Vec<usize>, SplitRejection); 3] = [
        (vec![cluster(4, 0)], vec![], SplitRejection::SingleCluster),
        (vec![], cluster(4, 0), SplitRejection::SingleCluster),
        (
            vec![cluster(2, 0), cluster(2, 2)],
            cluster(7, 4),
            SplitRejection::AllNoise,
        ),
    ];

    for (clusters, noise, expected) in cases {
        assert_eq!(judge_split(clusters, noise, size_of).unwrap_err(), expected);
    }
}

/// The noise cap is the port's own, not flight's, and it is the one that overrules a

/// A piece under the output floor is still a piece. Pouring it in with the noise denies it
/// the recruitment and merge passes that could carry it over the floor.
#[test]
fn a_piece_too_small_to_write_is_still_kept() {
    let (kept, spare) = judge_split(
        vec![cluster(3, 0), cluster(3, 3), cluster(1, 6)],
        cluster(1, 7),
        size_of,
    )
    .unwrap();

    assert_eq!(kept.len(), 3);
    assert_eq!(spare, vec![7]);
}

/// Noise is a piece too: one cluster beside it is a bin that lost something, not a bin left
/// whole.
#[test]
fn a_lone_cluster_beside_noise_is_a_split() {
    let (kept, spare) = judge_split(vec![cluster(3, 0)], cluster(1, 3), size_of).unwrap();

    assert_eq!(kept.len(), 1);
    assert_eq!(spare, vec![3]);
}

/// A bin that comes apart into one bin and dust has not been split, it has lost contigs,
/// and eject is the pass with a bar for that.
#[test]
fn one_survivor_beside_dust_is_a_shred() {
    let floor = 3 * CONTIG_LENGTH;
    assert!(!leaves_two_standing(
        &[cluster(5, 0), cluster(1, 5), cluster(1, 6)],
        floor,
        size_of
    ));
    assert!(leaves_two_standing(
        &[cluster(3, 0), cluster(1, 3), cluster(4, 4)],
        floor,
        size_of
    ));
}

const TIGHT: usize = 40;
const STRAY: usize = TIGHT;
const TIGHT_LENGTH: usize = 10_000;
const CLOSED_LENGTH: usize = 3_000_000;
const MIN_BIN_SIZE: usize = 200_000;

/// One bin of co-abundant, compositionally alike short contigs holding a long contig that
/// matches neither, which is the shape every absorbed closed genome takes.
fn absorbed_fixture() -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let n = TIGHT + 1;
    let mut coverage = Array2::zeros((n, 4));
    let mut tnf = Array2::zeros((n, 6));
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    let mut noise = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };

    for row in 0..TIGHT {
        coverage[[row, 0]] = 10.0 + noise();
        coverage[[row, 1]] = 2.0;
        coverage[[row, 2]] = 4.0 + noise();
        coverage[[row, 3]] = 1.0;
        for (column, base) in [0.1, -0.2, 0.3, -0.4, 0.2, -0.1].iter().enumerate() {
            tnf[[row, column]] = base + 0.1 * noise();
        }
    }
    coverage[[STRAY, 0]] = 60.0;
    coverage[[STRAY, 1]] = 5.0;
    coverage[[STRAY, 2]] = 0.5;
    coverage[[STRAY, 3]] = 0.2;
    for (column, base) in [-1.5, 1.2, -0.9, 1.4, -1.1, 0.8].iter().enumerate() {
        tnf[[STRAY, column]] = *base;
    }

    let mut lengths = vec![TIGHT_LENGTH; n];
    lengths[STRAY] = CLOSED_LENGTH;
    (coverage, tnf, lengths)
}

#[test]
fn a_closed_contig_comes_out_of_the_bin_that_absorbed_it() {
    let (coverage, tnf, lengths) = absorbed_fixture();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let objective = ObjectiveChoice::parse("codelength").unwrap().build();
    let settings = RefineSettings {
        min_bin_size: MIN_BIN_SIZE,
        max_bin_size: 15_000_000,
        n_neighbours: 100,
        max_retries: 5,
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
    let bins = BTreeMap::from([(0usize, (0..=TIGHT).collect::<Vec<_>>())]);
    let mut refiner = Refiner::new(features, &objective, settings, bins, Vec::new());

    assert!(refiner.run() >= 1);
    assert!(
        refiner.bins.values().any(|bin| bin == &vec![STRAY]),
        "{:?}",
        refiner.bins
    );
    assert!(refiner.unbinned.is_empty());
}
