use std::collections::HashMap;

fn pairs(count: f64) -> f64 {
    count * (count - 1.0) / 2.0
}

/// Chance corrected agreement between two labellings. Two identical trivial partitions still
/// score 1.0, so coarse resolutions are stable for free and this cannot rank them alone.
pub fn adjusted_rand_index(left: &[i32], right: &[i32]) -> f64 {
    if left.is_empty() || left.len() != right.len() {
        return f64::NAN;
    }

    let mut joint: HashMap<(i32, i32), f64> = HashMap::new();
    let mut left_totals: HashMap<i32, f64> = HashMap::new();
    let mut right_totals: HashMap<i32, f64> = HashMap::new();
    for (a, b) in left.iter().zip(right) {
        *joint.entry((*a, *b)).or_default() += 1.0;
        *left_totals.entry(*a).or_default() += 1.0;
        *right_totals.entry(*b).or_default() += 1.0;
    }

    let total = pairs(left.len() as f64);
    if total <= 0.0 {
        return f64::NAN;
    }

    let agreements = joint.values().copied().map(pairs).sum::<f64>();
    let left_pairs = left_totals.values().copied().map(pairs).sum::<f64>();
    let right_pairs = right_totals.values().copied().map(pairs).sum::<f64>();

    let expected = left_pairs * right_pairs / total;
    let largest = (left_pairs + right_pairs) / 2.0;
    if (largest - expected).abs() <= 1e-12 * largest.max(1.0) {
        return 1.0;
    }
    (agreements - expected) / (largest - expected)
}

pub fn mean_pairwise(labellings: &[Vec<i32>]) -> f64 {
    let mut total = 0.0;
    let mut compared = 0.0;
    for (index, left) in labellings.iter().enumerate() {
        for right in &labellings[index + 1..] {
            total += adjusted_rand_index(left, right);
            compared += 1.0;
        }
    }
    if compared == 0.0 {
        return f64::NAN;
    }
    total / compared
}
