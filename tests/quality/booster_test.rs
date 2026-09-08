//! Reading a boosted tree back out of its text form and walking it.

use rosella::quality::booster::Booster;

const MODEL: &str = "tree
version=v3
num_class=1
objective=regression sqrt


Tree=0
num_leaves=3
split_feature=0 2
threshold=1.5 10.5
decision_type=2 2
left_child=1 -1
right_child=-2 -3
leaf_value=0.5 2 1
shrinkage=1


Tree=1
num_leaves=1
leaf_value=0.25
shrinkage=0.1

end of trees
";

/// A child index is a node when it is positive and the ones complement of a leaf when it is
/// negative, which is the one thing a hand written walk gets wrong.
#[test]
fn a_negative_child_names_a_leaf() {
    let booster = Booster::parse(MODEL).expect("model");

    // feature 0 above the split goes right, to leaf 1.
    let right = booster.predict(&[2.0, 0.0, 0.0]);
    assert!((right - (2.0 + 0.25f64).powi(2)).abs() < 1e-12, "{right}");

    // below it, feature 2 decides between leaves 0 and 2.
    let low = booster.predict(&[1.0, 0.0, 5.0]);
    assert!((low - (0.5 + 0.25f64).powi(2)).abs() < 1e-12, "{low}");
    let high = booster.predict(&[1.0, 0.0, 50.0]);
    assert!((high - (1.0 + 0.25f64).powi(2)).abs() < 1e-12, "{high}");
}

/// The models were trained on the square root of the answer, so a negative sum has to come
/// back negative rather than squaring into a large positive completeness.
#[test]
fn a_negative_sum_stays_negative() {
    let booster = Booster::parse(
        "tree
Tree=0
num_leaves=1
leaf_value=-3
shrinkage=1

end of trees
",
    )
    .expect("model");
    assert!((booster.predict(&[0.0]) + 9.0).abs() < 1e-12);
}

#[test]
fn a_model_with_no_trees_is_refused() {
    assert!(Booster::parse("tree\nversion=v3\n").is_err());
}
