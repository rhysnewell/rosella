use std::collections::HashMap;
use std::hash::Hash;

use rand::Rng;

pub fn positions<K: Hash + Eq>(keys: impl IntoIterator<Item = K>) -> HashMap<K, usize> {
    keys.into_iter()
        .enumerate()
        .map(|(at, key)| (key, at))
        .collect()
}

pub fn draw<'a, T>(rng: &mut impl Rng, items: &'a mut [T], at: usize) -> &'a T {
    let pick = rng.random_range(at..items.len());
    items.swap(at, pick);
    &items[at]
}

pub fn reorder<T: Clone>(values: &[T], order: &[usize]) -> Vec<T> {
    order.iter().map(|at| values[*at].clone()).collect()
}
