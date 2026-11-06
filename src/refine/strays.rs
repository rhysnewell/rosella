use std::cmp::Reverse;
use std::collections::{BTreeMap, HashMap, HashSet};

use crate::quality::Scorer;
use crate::refine::rung::{Bars, RUNGS};

fn find(parent: &mut HashMap<usize, usize>, contig: usize) -> usize {
    let mut at = contig;
    loop {
        let up = *parent.entry(at).or_insert(at);
        if up == at {
            return at;
        }
        let above = *parent.entry(up).or_insert(up);
        parent.insert(at, above);
        at = above;
    }
}

// Each chain's root is its smallest contig, so the chains come out in one order on every run.
fn chains(bridges: &[(usize, usize, u32)]) -> Vec<Vec<usize>> {
    let mut parent = HashMap::new();
    for (one, other, _) in bridges {
        let (one, other) = (find(&mut parent, *one), find(&mut parent, *other));
        parent.insert(one.max(other), one.min(other));
    }
    let mut chains = BTreeMap::<usize, Vec<usize>>::new();
    for contig in parent.keys().copied().collect::<Vec<_>>() {
        chains
            .entry(find(&mut parent, contig))
            .or_default()
            .push(contig);
    }
    chains
        .into_values()
        .map(|mut chain| {
            chain.sort_unstable();
            chain
        })
        .collect()
}

// A bridge between two bins mostly touches one stray contig, so moving the strays home keeps the
// rest of the other bin where it is. Merging the whole bins drags a second genome in.
pub fn gather(
    bins: &mut BTreeMap<usize, Vec<usize>>,
    unbinned: &mut Vec<usize>,
    bridges: &[(usize, usize, u32)],
    lengths: &[usize],
    quality: &dyn Scorer,
    bars: &Bars,
) -> Vec<(usize, Option<usize>, usize)> {
    let mut owner = bins
        .iter()
        .flat_map(|(bin, contigs)| contigs.iter().map(move |contig| (*contig, Some(*bin))))
        .chain(unbinned.iter().map(|contig| (*contig, None)))
        .collect::<HashMap<_, _>>();
    let held = bridges
        .iter()
        .filter(|(one, other, _)| owner.contains_key(one) && owner.contains_key(other))
        .copied()
        .collect::<Vec<_>>();
    let mut moved = Vec::new();
    for chain in chains(&held) {
        let mut bases = BTreeMap::<usize, usize>::new();
        for contig in &chain {
            if let Some(bin) = owner[contig] {
                *bases.entry(bin).or_default() += lengths[*contig];
            }
        }
        let Some(home) = bases
            .into_iter()
            .max_by_key(|(bin, held)| (*held, Reverse(*bin)))
            .map(|(bin, _)| bin)
        else {
            continue;
        };
        let mut strays = chain
            .into_iter()
            .filter(|contig| owner[contig] != Some(home))
            .collect::<Vec<_>>();
        let held = anchored(&strays, &owner, bins, quality, bars);
        strays.retain(|contig| !held.contains(contig));
        if strays.is_empty() {
            continue;
        }
        let members = &bins[&home];
        let mut grown = members.iter().chain(&strays).copied().collect::<Vec<_>>();
        if quality.score(&grown).contamination > quality.score(members).contamination {
            continue;
        }
        for stray in strays {
            let from = owner.insert(stray, Some(home)).flatten();
            if let Some(members) = from.and_then(|bin| bins.get_mut(&bin)) {
                members.retain(|contig| *contig != stray);
            }
            moved.push((stray, from, home));
        }
        grown.sort_unstable();
        bins.insert(home, grown);
    }
    unbinned.retain(|contig| owner[contig].is_none());
    bins.retain(|_, contigs| !contigs.is_empty());
    moved
}

// The marker gate on the home bin cannot see what a stray's own bin loses. A bin on the rung
// ladder is a genome already, so it keeps every stray its completeness rests on.
fn anchored(
    strays: &[usize],
    owner: &HashMap<usize, Option<usize>>,
    bins: &BTreeMap<usize, Vec<usize>>,
    quality: &dyn Scorer,
    bars: &Bars,
) -> HashSet<usize> {
    let mut leaving = BTreeMap::<usize, HashSet<usize>>::new();
    for stray in strays {
        if let Some(bin) = owner[stray] {
            leaving.entry(bin).or_default().insert(*stray);
        }
    }
    leaving
        .into_iter()
        .filter(|(bin, out)| {
            let members = &bins[bin];
            let left = members
                .iter()
                .copied()
                .filter(|contig| !out.contains(contig))
                .collect::<Vec<_>>();
            let before = quality.score(members);
            bars.rung_of(before) < RUNGS && quality.score(&left).completeness < before.completeness
        })
        .flat_map(|(_, out)| out)
        .collect()
}
