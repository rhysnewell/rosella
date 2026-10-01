use std::cmp::Ordering;
use std::collections::HashSet;

// Worth decides, and the smaller candidate breaks a tie so a heap of equally worthy
// proposals drains in one order rather than in hash order. `extra` is payload, never compared.
pub struct Ranked<T> {
    pub worth: f64,
    pub contigs: Vec<usize>,
    pub extra: T,
}

impl<T> PartialEq for Ranked<T> {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl<T> Eq for Ranked<T> {}

impl<T> PartialOrd for Ranked<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T> Ord for Ranked<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.worth
            .total_cmp(&other.worth)
            .then_with(|| other.contigs.cmp(&self.contigs))
    }
}

// A proposal the winners already emptied is not the proposal that was scored, so what is left
// of it is put back through the same bar rather than trusted on the rank it earned whole.
pub fn remaining(contigs: &[usize], claimed: &HashSet<usize>) -> Vec<usize> {
    contigs
        .iter()
        .copied()
        .filter(|contig| !claimed.contains(contig))
        .collect()
}

pub fn remaining_in(contigs: &[usize], pool: &HashSet<usize>) -> Vec<usize> {
    contigs
        .iter()
        .copied()
        .filter(|contig| pool.contains(contig))
        .collect()
}

// Contig order decides bin ids downstream, so every set that becomes a bin is sorted here
// rather than left in hash order.
pub fn sorted(contigs: impl IntoIterator<Item = usize>) -> Vec<usize> {
    let mut contigs = contigs.into_iter().collect::<Vec<_>>();
    contigs.sort_unstable();
    contigs
}
