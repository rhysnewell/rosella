use std::collections::{BTreeMap, HashMap};

use rayon::prelude::*;

use crate::kmers::sketch::ContigSketches;

// A genome's pieces share no sequence, so a contig its bin-mates mostly hold again is one more
// copy of an element the assembly kept several times, not a piece of the genome.
const COPY_SHARE: f64 = 0.5;

pub fn copies(bins: &BTreeMap<usize, Vec<usize>>, sketches: &ContigSketches) -> Vec<Vec<usize>> {
    bins.par_iter()
        .flat_map_iter(|(_, members)| elements(members, sketches))
        .collect()
}

// A repeat inside a longer genome contig is held by a contig that is no copy itself, so only
// copies whose best match is also a copy leave, and they leave together.
fn elements(members: &[usize], sketches: &ContigSketches) -> Vec<Vec<usize>> {
    let mut holders = HashMap::<u64, Vec<usize>>::new();
    for (at, contig) in members.iter().enumerate() {
        for hash in sketches.sketch(*contig) {
            holders.entry(*hash).or_default().push(at);
        }
    }
    let partner = members
        .iter()
        .enumerate()
        .map(|(at, contig)| best(at, sketches.sketch(*contig), &holders))
        .collect::<Vec<_>>();
    let mut root = (0..members.len()).collect::<Vec<_>>();
    for (at, found) in partner.iter().enumerate() {
        if let Some(other) = found
            && partner[*other].is_some()
        {
            let (left, right) = (find(&mut root, at), find(&mut root, *other));
            root[left.max(right)] = left.min(right);
        }
    }
    let mut groups = BTreeMap::<usize, Vec<usize>>::new();
    for (at, contig) in members.iter().enumerate() {
        if partner[at].is_some_and(|other| partner[other].is_some()) {
            groups.entry(find(&mut root, at)).or_default().push(*contig);
        }
    }
    groups
        .into_values()
        .filter(|group| group.len() > 1)
        .collect()
}

fn best(at: usize, hashes: &[u64], holders: &HashMap<u64, Vec<usize>>) -> Option<usize> {
    let mut shared = 0;
    let mut with = HashMap::<usize, usize>::new();
    for hash in hashes {
        let others = holders[hash].iter().filter(|other| **other != at);
        let mut held = false;
        for other in others {
            *with.entry(*other).or_default() += 1;
            held = true;
        }
        shared += usize::from(held);
    }
    (!hashes.is_empty() && shared as f64 > COPY_SHARE * hashes.len() as f64)
        .then(|| {
            with.into_iter()
                .max_by_key(|(other, count)| (*count, usize::MAX - other))
        })
        .flatten()
        .map(|(other, _)| other)
}

fn find(root: &mut [usize], at: usize) -> usize {
    let mut top = at;
    while root[top] != top {
        top = root[top];
    }
    let mut walk = at;
    while root[walk] != top {
        let next = root[walk];
        root[walk] = top;
        walk = next;
    }
    top
}
