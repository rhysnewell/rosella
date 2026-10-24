use std::collections::HashMap;

use crate::recover::origin::BinSet;
use crate::refine::owners::heir;

#[derive(Default, Debug, PartialEq)]
pub struct Ballot {
    pub bins: Vec<usize>,
    pub voters: usize,
}

// An offered bin speaks for the final bin holding most of its long bases, since long contigs
// decide bins and short ones follow. A set that bins nothing as short as a contig has no say.
pub fn ballots(
    order: &[usize],
    lengths: &[usize],
    cutoff: usize,
    sets: &[BinSet],
    bin_of: &HashMap<usize, usize>,
) -> Vec<Ballot> {
    if sets.is_empty() {
        return Vec::new();
    }
    let slot = order
        .iter()
        .enumerate()
        .map(|(at, contig)| (*contig, at))
        .collect::<HashMap<_, _>>();
    let mut ballots = order.iter().map(|_| Ballot::default()).collect::<Vec<_>>();
    let mut shortest = Vec::with_capacity(sets.len());
    for set in sets {
        let members = set.bins.iter().flat_map(|(_, members)| members);
        shortest.extend(members.map(|contig| lengths[*contig]).min());
        for (_, members) in &set.bins {
            let long = members
                .iter()
                .copied()
                .filter(|contig| lengths[*contig] >= cutoff)
                .collect::<Vec<_>>();
            let Some((bin, held)) = heir(&long, bin_of, |contig| lengths[contig]) else {
                continue;
            };
            if 2 * held <= long.iter().map(|contig| lengths[*contig]).sum::<usize>() {
                continue;
            }
            for at in members.iter().filter_map(|contig| slot.get(contig)) {
                ballots[*at].bins.push(bin);
            }
        }
    }
    shortest.sort_unstable();
    for (contig, ballot) in order.iter().zip(&mut ballots) {
        ballot.voters = shortest.partition_point(|least| *least <= lengths[*contig]);
    }
    ballots
}
