//! The sweep ranks labellings with this, and the writer then keeps them on assembled bp.
//! These cover the two knobs that bring the score into the writer's unit, and the reason
//! the objective samples for itself rather than being handed a sample.

use ndarray::Array2;
use rosella::clustering::objective::{
    ClusterObjective, ClusterWeight, Dbcv, EmbeddingSample, ObjectiveChoice,
};

const SEED: u64 = 42;
const SEPARATION: f64 = 100.0;

/// Two well separated blobs, laid out without an RNG so the geometry is fixed.
fn embedding(sizes: &[usize]) -> (Array2<f64>, Vec<i32>) {
    let total = sizes.iter().sum::<usize>();
    let mut points = Array2::zeros((total, 2));
    let mut labels = Vec::with_capacity(total);
    let mut row = 0;
    for (cluster, size) in sizes.iter().enumerate() {
        for index in 0..*size {
            points[[row, 0]] = cluster as f64 * SEPARATION + (index % 40) as f64 * 0.03;
            points[[row, 1]] = (index / 40) as f64 * 0.03;
            labels.push(cluster as i32);
            row += 1;
        }
    }
    (points, labels)
}

fn contigs(n: usize) -> Vec<usize> {
    (0..n).collect()
}

fn score(lengths: &[usize], weight: ClusterWeight, floor: Option<usize>, sizes: &[usize]) -> f64 {
    let (points, labels) = embedding(sizes);
    let sample = EmbeddingSample::new(points.view(), SEED);
    Dbcv::new(lengths)
        .with_weight(weight)
        .with_floor(floor)
        .score(&sample, &contigs(labels.len()), &labels)
}

/// Cluster 0 is many short contigs, cluster 1 is few long ones. Counting contigs makes the
/// first almost the whole score and counting bases makes the second.
fn lopsided() -> (Vec<usize>, Vec<usize>) {
    let sizes = vec![120, 30];
    let mut lengths = vec![1_000; 120];
    lengths.extend(std::iter::repeat_n(100_000, 30));
    (sizes, lengths)
}

#[test]
fn bp_and_count_weighting_disagree_when_bp_share_is_not_contig_share() {
    let (sizes, lengths) = lopsided();
    let counted = score(&lengths, ClusterWeight::Count, None, &sizes);
    let weighed = score(&lengths, ClusterWeight::Bp, None, &sizes);

    assert!(counted > 0.0 && weighed > 0.0);
    assert!(
        (counted - weighed).abs() > 1e-6,
        "count {counted} and bp {weighed} agreed, so the weight is not reaching the score"
    );
}

#[test]
fn the_default_ignores_the_lengths_it_is_handed() {
    let (sizes, lengths) = lopsided();
    let uniform = vec![1_000; lengths.len()];
    assert_eq!(
        score(&lengths, ClusterWeight::Count, None, &sizes),
        score(&uniform, ClusterWeight::Count, None, &sizes),
    );
}

#[test]
fn the_floor_drops_a_cluster_under_it_and_leaves_one_over_it() {
    let (sizes, lengths) = lopsided();
    // Cluster 0 is 120 kb, cluster 1 is 3 Mb.
    let unfloored = score(&lengths, ClusterWeight::Count, None, &sizes);
    let over_the_small = score(&lengths, ClusterWeight::Count, Some(200_000), &sizes);
    let over_both = score(&lengths, ClusterWeight::Count, Some(4_000_000), &sizes);
    let under_both = score(&lengths, ClusterWeight::Count, Some(100_000), &sizes);

    assert_eq!(
        under_both, unfloored,
        "a floor below every cluster is inert"
    );
    assert!(over_the_small < unfloored && over_the_small > 0.0);
    assert_eq!(over_both, 0.0, "a floor above every cluster leaves nothing");
}

/// The objective is handed every row and samples for itself. Were it handed a sample, a
/// cluster's bp would be the sample's bp and a floor test would reject bins that are whole.
#[test]
fn bp_is_read_off_the_whole_labelling_rather_than_the_sample() {
    let sizes = vec![5_000, 5_000];
    let lengths = vec![100_000; 10_000];
    // Each cluster is 500 Mb whole and about 250 Mb once halved by the sample limit.
    let floor = Some(400_000_000);

    let unfloored = score(&lengths, ClusterWeight::Count, None, &sizes);
    assert!(
        unfloored > 0.0,
        "nothing scored, so the comparison below is vacuous"
    );
    assert_eq!(
        score(&lengths, ClusterWeight::Count, floor, &sizes),
        unfloored,
        "the floor rejected a cluster that is over it, so it was applied to the sample"
    );
}

#[test]
fn every_objective_name_resolves_and_an_unknown_one_does_not() {
    for (name, weight, floor) in [
        ("dbcv", ClusterWeight::Count, false),
        ("dbcv-bp", ClusterWeight::Bp, false),
        ("dbcv-floor", ClusterWeight::Count, true),
        ("dbcv-bp-floor", ClusterWeight::Bp, true),
    ] {
        let choice = ObjectiveChoice::parse(name).expect(name);
        assert_eq!(choice.weight, weight);
        assert_eq!(choice.floor, floor);
    }
    assert!(ObjectiveChoice::parse("markers").is_none());
}
