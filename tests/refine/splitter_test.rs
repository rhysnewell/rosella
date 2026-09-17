//! Whether the re-clustering that came back is worth taking. Whether the bin was worth
//! re-clustering in the first place is `bar_test`.

use std::collections::BTreeMap;

use ndarray::Array2;
use rosella::clustering::graph_partition::Partition;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::gates::SplitRejection;
use rosella::refine::proposal::{Standing, judge_split, standing};
use rosella::refine::splitter::{RefineSettings, Refiner};
use rosella::seeds::Seeds;

#[path = "../support/scorer.rs"]
mod scorer;

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

/// A piece under the output floor is still a piece. Pouring it in with the noise denies it
/// the merge pass that could carry it over the floor.
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

#[test]
fn one_survivor_beside_dust_is_a_shred() {
    let floor = 3 * CONTIG_LENGTH;
    assert!(!matches!(
        standing(
            &[cluster(5, 0), cluster(1, 5), cluster(1, 6)],
            0,
            floor,
            size_of
        ),
        Standing::Many
    ));
    assert!(matches!(
        standing(
            &[cluster(3, 0), cluster(1, 3), cluster(4, 4)],
            0,
            floor,
            size_of
        ),
        Standing::Many
    ));
}

/// Trimming takes the same cut, so the bar moves to what walks away. A piece too small to be
/// written as a bin is dust; anything that could have been a bin is a genome coming apart.
#[test]
fn a_trim_turns_on_what_leaves_not_what_stands() {
    let floor = 3 * CONTIG_LENGTH;
    let bin_floor = 2 * CONTIG_LENGTH;
    let walks = |pieces: &[Vec<usize>], scattered| matches!(standing(pieces, scattered, floor, size_of), Standing::One { largest } if largest < bin_floor);
    assert!(walks(&[cluster(5, 0), cluster(1, 5), cluster(1, 6)], 0));
    assert!(!walks(&[cluster(5, 0), cluster(2, 5), cluster(2, 7)], 0));
    assert!(!walks(&[cluster(5, 0), cluster(1, 5)], bin_floor));
}

const TIGHT: usize = 40;
const STRAY: usize = TIGHT;
const TIGHT_LENGTH: usize = 10_000;
const CLOSED_LENGTH: usize = 3_000_000;
const MIN_BIN_SIZE: usize = 200_000;

/// One bin of co-abundant, compositionally alike short contigs holding a long contig that
/// matches neither, which is the shape every absorbed closed genome takes. `anchor` adds a
/// second long contig for its own bin, which is what gives the run a measured genome floor.
fn absorbed_fixture(anchor: Option<usize>, tight: usize) -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let n = TIGHT + 1 + usize::from(anchor.is_some());
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

    let mut lengths = vec![tight; n];
    lengths[STRAY] = CLOSED_LENGTH;
    if let Some(length) = anchor {
        let row = n - 1;
        coverage[[row, 0]] = 120.0;
        coverage[[row, 1]] = 9.0;
        coverage[[row, 2]] = 0.1;
        coverage[[row, 3]] = 0.05;
        for (column, base) in [1.4, -1.3, 1.1, -1.2, 0.9, -0.7].iter().enumerate() {
            tnf[[row, column]] = *base;
        }
        lengths[row] = length;
    }
    (coverage, tnf, lengths)
}

fn settings(trim: bool) -> RefineSettings {
    RefineSettings {
        min_bin_size: MIN_BIN_SIZE,
        max_bin_size: 15_000_000,
        knn_candidates: rosella::embedding::knn::MAX_CANDIDATES,
        n_neighbours: 100,
        max_retries: 5,
        seeds: Seeds {
            knn: 42,
            seed: 42,
            partition: 42,
        },
        max_contamination: None,
        partition: Partition::Both,
        trim,
        anchor_ladder: false,
        leiden: rosella::clustering::leiden::Null::default(),
    }
}

#[test]
fn a_closed_contig_comes_out_of_the_bin_that_absorbed_it() {
    let (coverage, tnf, lengths) = absorbed_fixture(None, TIGHT_LENGTH);
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let bins = BTreeMap::from([(0usize, (0..=TIGHT).collect::<Vec<_>>())]);
    let mut refiner = Refiner::new(features, settings(false), bins, Vec::new());

    assert!(refiner.run() >= 1);
    assert!(
        refiner.bins.values().any(|bin| bin == &vec![STRAY]),
        "{:?}",
        refiner.bins
    );
    assert!(refiner.unbinned.is_empty());
}

/// The peel took an accepted outcome with no gate at all, so it could cut a genome off the
/// markers that left with the contig. Families shared with the peeled contig mean two
/// organisms; dust has no families to judge it by and is the peel working.
#[test]
fn a_peel_is_refused_when_it_cuts_a_genome_off_its_markers() {
    let dust = MIN_BIN_SIZE / (2 * TIGHT);
    let cases = [
        (TIGHT_LENGTH, Some(7), Some(7), true),
        (TIGHT_LENGTH, Some(7), Some(9), false),
        (TIGHT_LENGTH, Some(7), None, false),
        (TIGHT_LENGTH, None, None, true),
        (dust, Some(7), None, true),
    ];

    for (tight, peeled, held, splits) in cases {
        let (coverage, tnf, lengths) = absorbed_fixture(None, tight);
        let mut families = vec![held; lengths.len()];
        families[STRAY] = peeled;
        let quality = scorer::FamilyScorer::new(families);
        let features = ContigFeatures::new(&coverage, &tnf, &lengths);
        let absorbed = (0..=TIGHT).collect::<Vec<_>>();
        let bins = BTreeMap::from([(0usize, absorbed.clone())]);
        let mut refiner =
            Refiner::new(features, settings(false), bins, Vec::new()).with_quality(&quality);
        refiner.run();

        assert_eq!(
            refiner.bins.values().any(|bin| bin == &vec![STRAY]),
            splits,
            "tight {tight}, peeled {peeled:?}, held {held:?}: {:?}",
            refiner.bins
        );
    }
}
