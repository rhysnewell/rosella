use std::collections::HashMap;
use std::ops::Range;

use crate::recover::origin::BinSet;
use crate::refine::owners::{credit, heir};

const FREE: u32 = u32::MAX;
const LOOSE: u32 = u32::MAX - 1;

#[derive(Default, Debug, PartialEq)]
pub struct Ballot {
    pub bins: Vec<(usize, f32)>,
    pub say: f32,
}

impl Ballot {
    // The neighbourhood is one voter, so a contig no bin can win a majority of never joins.
    pub fn can_join(&self) -> bool {
        let mut tally = Vec::new();
        for (bin, weight) in &self.bins {
            credit(&mut tally, *bin, *weight);
        }
        let most = tally.iter().map(|(_, weight)| *weight).fold(0.0, f32::max);
        2.0 * (1.0 + most) > 1.0 + self.say
    }
}

// An offered bin speaks for the final bin holding most of its long bases, since long contigs
// decide bins and short ones follow. A set's say in a band is the share of the band it bins, from
// its own shortest contig up, so a set that rarely bins contigs that short rarely counts.
pub fn ballots(
    order: &[usize],
    spans: &[Range<usize>],
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
    for set in sets {
        let mut home = vec![FREE; order.len()];
        let mut shortest = usize::MAX;
        for (_, members) in &set.bins {
            shortest = members
                .iter()
                .map(|contig| lengths[*contig])
                .fold(shortest, usize::min);
            let target = majority(members, lengths, cutoff, bin_of).map_or(LOOSE, |bin| bin as u32);
            for at in members.iter().filter_map(|contig| slot.get(contig)) {
                home[*at] = target;
            }
        }
        for span in spans {
            let reach = span
                .clone()
                .filter(|at| lengths[order[*at]] >= shortest)
                .collect::<Vec<_>>();
            let binned = reach.iter().filter(|at| home[**at] != FREE).count();
            let rate = binned as f32 / reach.len().max(1) as f32;
            for at in reach {
                ballots[at].say += rate;
                if home[at] < LOOSE {
                    ballots[at].bins.push((home[at] as usize, rate));
                }
            }
        }
    }
    ballots
}

fn majority(
    members: &[usize],
    lengths: &[usize],
    cutoff: usize,
    bin_of: &HashMap<usize, usize>,
) -> Option<usize> {
    let long = members
        .iter()
        .copied()
        .filter(|contig| lengths[*contig] >= cutoff)
        .collect::<Vec<_>>();
    let (bin, held) = heir(&long, bin_of, |contig| lengths[contig])?;
    let total = long.iter().map(|contig| lengths[*contig]).sum::<usize>();
    (2 * held > total).then_some(bin)
}
