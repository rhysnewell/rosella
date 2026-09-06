use std::collections::BTreeMap;

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::refine::merger::{MergeBar, MergeSettings, merge_bins};

const CONTIG_LENGTH: usize = 100_000;
const PER_BIN: usize = 6;
const GENOME_LENGTH: usize = 1_500_000;

fn settings(genome_floor: Option<usize>) -> MergeSettings {
    MergeSettings {
        genome_floor,
        max_bin_size: usize::MAX,
        seed: 42,
        ..MergeSettings::default()
    }
}

fn jitter(seed: u64) -> impl FnMut() -> f64 {
    let mut state = seed;
    move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn two_clouds() -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let rows = PER_BIN * 2;
    let mut coverage = Array2::zeros((rows, 2));
    let mut tnf = Array2::zeros((rows, 4));
    let mut next = jitter(0x2545_F491_4F6C_DD1D);

    for row in 0..rows {
        let second = if row < PER_BIN { 0.0 } else { 1.0 };
        coverage[[row, 0]] = 10.0 + second * 4.0 + 6.0 * next();
        coverage[[row, 1]] = 4.0 + next();
        for column in 0..4 {
            tnf[[row, column]] = 0.1 * column as f64 + second * 0.2 + 0.4 * next();
        }
    }
    (coverage, tnf, vec![CONTIG_LENGTH; rows])
}

/// Centroids 0.308 apart against a bar of 0.345, so the first gate accepts, and a union at
/// 0.395 that the second refuses. Removing that refusal costs CAMI I high 45 bins at t5, so
/// the rejection is load-bearing however badly the two quantities match.
#[test]
fn the_second_gate_refuses_what_the_first_accepts() {
    let (coverage, tnf, lengths) = two_clouds();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let (merged, merges) = merge_bins(
        &features,
        BTreeMap::from([
            (0, (0..PER_BIN).collect::<Vec<_>>()),
            (1, (PER_BIN..PER_BIN * 2).collect()),
        ]),
        settings(None),
    );

    assert_eq!(merges, 0);
    assert_eq!(merged.len(), 2);
}

fn one_cloud_and_a_singleton() -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let rows = PER_BIN * 2 + 1;
    let mut coverage = Array2::zeros((rows, 2));
    let mut tnf = Array2::zeros((rows, 4));
    let mut next = jitter(0x2545_F491_4F6C_DD1D);

    for row in 0..rows {
        coverage[[row, 0]] = 10.0 + 0.5 * next();
        coverage[[row, 1]] = 4.0 + 0.5 * next();
        for column in 0..4 {
            tnf[[row, column]] = 0.1 * column as f64 + 0.05 * next();
        }
    }
    (coverage, tnf, vec![CONTIG_LENGTH; rows])
}

/// Without a genome floor a bin of one contig has no spread to be judged against, so it is
/// carried through the pass rather than dropped. This is what `--no-merge-singles` leaves.
#[test]
fn a_one_contig_bin_survives_a_merge_it_cannot_join() {
    let (coverage, tnf, lengths) = one_cloud_and_a_singleton();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let (merged, merges) = merge_bins(
        &features,
        BTreeMap::from([
            (0, (0..PER_BIN).collect::<Vec<_>>()),
            (1, (PER_BIN..PER_BIN * 2).collect()),
            (2, vec![PER_BIN * 2]),
        ]),
        settings(None),
    );

    assert_eq!(merges, 1);
    assert_eq!(merged.len(), 2);
    assert_eq!(merged[&2], vec![PER_BIN * 2]);
    assert_eq!(merged[&0].len(), PER_BIN * 2);
}

/// The contig sits in the same cloud as both bins, so the spread they already tolerate holds
/// it and the three collapse into one.
#[test]
fn a_one_contig_bin_joins_the_bin_whose_spread_holds_it() {
    let (coverage, tnf, lengths) = one_cloud_and_a_singleton();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let (merged, merges) = merge_bins(
        &features,
        BTreeMap::from([
            (0, (0..PER_BIN).collect::<Vec<_>>()),
            (1, (PER_BIN..PER_BIN * 2).collect()),
            (2, vec![PER_BIN * 2]),
        ]),
        settings(Some(GENOME_LENGTH)),
    );

    assert_eq!(merges, 2);
    assert_eq!(merged.len(), 1);
    assert_eq!(merged[&0].len(), PER_BIN * 2 + 1);
}

#[test]
fn two_one_contig_bins_have_no_spread_to_judge() {
    let (coverage, tnf, lengths) = one_cloud_and_a_singleton();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let (merged, merges) = merge_bins(
        &features,
        BTreeMap::from([(0, vec![0]), (1, vec![1])]),
        settings(Some(GENOME_LENGTH)),
    );

    assert_eq!(merges, 0);
    assert_eq!(merged.len(), 2);
}

fn a_cloud_with_two_long_contigs() -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let rows = PER_BIN + 1;
    let mut coverage = Array2::zeros((rows, 2));
    let mut tnf = Array2::zeros((rows, 4));
    let mut next = jitter(0x2545_F491_4F6C_DD1D);

    for row in 0..rows {
        coverage[[row, 0]] = 10.0 + 0.5 * next();
        coverage[[row, 1]] = 4.0 + 0.5 * next();
        for column in 0..4 {
            tnf[[row, column]] = 0.1 * column as f64 + 0.05 * next();
        }
    }
    let mut lengths = vec![CONTIG_LENGTH; rows];
    lengths[0] = GENOME_LENGTH;
    lengths[PER_BIN] = GENOME_LENGTH;
    (coverage, tnf, lengths)
}

/// Solo stands two genome-sized contigs apart, so a merge that puts them back together undoes
/// it. The same input joins once the floor that says what genome-sized means is raised past it.
#[test]
fn a_merge_that_solo_would_take_apart_is_refused() {
    let (coverage, tnf, lengths) = a_cloud_with_two_long_contigs();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let bins = BTreeMap::from([(0, (0..PER_BIN).collect::<Vec<_>>()), (1, vec![PER_BIN])]);

    let (_, vetoed) = merge_bins(&features, bins.clone(), settings(Some(1_000_000)));
    let (_, allowed) = merge_bins(&features, bins, settings(Some(GENOME_LENGTH * 4)));

    assert_eq!(vetoed, 0);
    assert_eq!(allowed, 1);
}

fn widest(genome_floor: Option<usize>) -> MergeSettings {
    MergeSettings {
        bar: MergeBar::Widest,
        ..settings(genome_floor)
    }
}

/// The pair bar is the mean spread of two tight bins, so it bars exactly the clean halves it
/// should join. The loosest contig each bin already holds is a scale it has agreed to.
#[test]
fn the_widest_bar_joins_what_the_pair_bar_refuses() {
    let (coverage, tnf, lengths) = two_clouds();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let bins = BTreeMap::from([
        (0, (0..PER_BIN).collect::<Vec<_>>()),
        (1, (PER_BIN..PER_BIN * 2).collect()),
    ]);

    let (_, refused) = merge_bins(&features, bins.clone(), settings(None));
    let (merged, joined) = merge_bins(&features, bins, widest(None));

    assert_eq!(refused, 0);
    assert_eq!(joined, 1);
    assert_eq!(merged.len(), 1);
}

#[test]
fn a_complete_bin_does_not_recruit() {
    let (coverage, tnf, lengths) = two_clouds();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let bins = BTreeMap::from([
        (0, (0..PER_BIN).collect::<Vec<_>>()),
        (1, (PER_BIN..PER_BIN * 2).collect()),
    ]);
    let short = |floor| MergeSettings {
        short_side: true,
        ..widest(Some(floor))
    };

    let (_, complete) = merge_bins(&features, bins.clone(), short(200_000));
    let (_, wanting) = merge_bins(&features, bins, short(500_000));

    assert_eq!(complete, 0);
    assert_eq!(wanting, 1);
}

fn three_clouds() -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let rows = PER_BIN * 3;
    let mut coverage = Array2::zeros((rows, 2));
    let mut tnf = Array2::zeros((rows, 4));
    let mut next = jitter(0x2545_F491_4F6C_DD1D);

    for row in 0..rows {
        let cloud = (row / PER_BIN) as f64;
        coverage[[row, 0]] = 10.0 + cloud * 2.0 + 4.0 * next();
        coverage[[row, 1]] = 4.0 + next();
        for column in 0..4 {
            tnf[[row, column]] = 0.1 * column as f64 + cloud * 0.05 + 0.3 * next();
        }
    }
    (coverage, tnf, vec![CONTIG_LENGTH; rows])
}

/// A bin whose nearest partner already has a nearer one of its own is not that partner's other
/// half, however close it sits. Reciprocity asks for that without naming a distance.
#[test]
fn only_a_pair_that_picks_each_other_merges() {
    let (coverage, tnf, lengths) = three_clouds();
    let features = ContigFeatures::new(&coverage, &tnf, &lengths);
    let bins = BTreeMap::from([
        (0, (0..PER_BIN).collect::<Vec<_>>()),
        (1, (PER_BIN..PER_BIN * 2).collect()),
        (2, (PER_BIN * 2..PER_BIN * 3).collect()),
    ]);

    let (_, all) = merge_bins(&features, bins.clone(), widest(None));
    let (merged, reciprocal) = merge_bins(
        &features,
        bins,
        MergeSettings {
            mutual: true,
            ..widest(None)
        },
    );

    assert_eq!(all, 2);
    assert_eq!(reciprocal, 1);
    assert_eq!(merged.len(), 2);
}
