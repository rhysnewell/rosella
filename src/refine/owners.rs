use std::cmp::Reverse;
use std::collections::HashMap;

pub fn owners<'a, B, C>(bins: B) -> HashMap<usize, usize>
where
    B: IntoIterator<Item = (usize, C)>,
    C: IntoIterator<Item = &'a usize>,
{
    bins.into_iter()
        .flat_map(|(label, members)| members.into_iter().map(move |contig| (*contig, label)))
        .collect()
}

// The lower label wins a tie, so the heir never depends on the order the bases were counted in.
pub fn heir(
    members: &[usize],
    owner: &HashMap<usize, usize>,
    length: impl Fn(usize) -> usize,
) -> Option<(usize, usize)> {
    let mut held = HashMap::<usize, usize>::new();
    for contig in members {
        if let Some(bin) = owner.get(contig) {
            *held.entry(*bin).or_default() += length(*contig);
        }
    }
    held.into_iter()
        .max_by_key(|(bin, bp)| (*bp, Reverse(*bin)))
}

// A neighbourhood names a handful of bins, so a short list beats a map built per contig.
pub fn credit<T: std::ops::AddAssign>(tally: &mut Vec<(usize, T)>, bin: usize, weight: T) {
    match tally.iter_mut().find(|(held, _)| *held == bin) {
        Some((_, held)) => *held += weight,
        None => tally.push((bin, weight)),
    }
}
