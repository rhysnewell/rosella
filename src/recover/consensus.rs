use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use rayon::prelude::*;

use crate::recover::ladder::{Judge, Offered};
use crate::refine::ranking::{remaining, sorted};

// Markers read a clean part of a genome as whole, so a winner's contigs are settled by which
// bins the other binners keep them in instead.
pub(crate) struct Ballots {
    voters: Vec<BTreeSet<usize>>,
    count: usize,
}

struct Rival {
    left: Vec<usize>,
    passes: bool,
}

impl Ballots {
    pub(crate) fn new(arms: &[Vec<usize>], voters: &[usize]) -> Self {
        Self {
            voters: arms
                .iter()
                .map(|held| held.iter().map(|arm| voters[*arm]).collect())
                .collect(),
            count: voters.iter().collect::<HashSet<_>>().len(),
        }
    }

    pub(crate) fn resolve(
        &self,
        judge: &Judge,
        offered: &Offered,
        at: usize,
        taken: &[usize],
        claimed: &HashSet<usize>,
    ) -> Option<Vec<usize>> {
        if self.count < 2 || !judge.passes(taken) {
            return None;
        }
        let primary = taken.iter().copied().collect::<HashSet<_>>();
        let rivals = taken
            .iter()
            .flat_map(|contig| offered.holders[contig].iter().copied())
            .filter(|rival| *rival != at)
            .collect::<BTreeSet<_>>()
            .into_par_iter()
            .map(|rival| {
                let left = remaining(&offered.found[rival], claimed);
                let passes = !left.is_empty() && judge.passes(&left);
                (rival, Rival { left, passes })
            })
            .collect::<HashMap<_, _>>();
        let own = &self.voters[at];
        let members = self.members(judge, taken, &primary, &rivals, own);
        if members.is_empty() {
            return None;
        }
        let size = own.len() + members.len();
        let held = members
            .iter()
            .map(|(voter, rival)| (*voter, rivals[rival].left.iter().copied().collect()))
            .collect::<BTreeMap<usize, HashSet<usize>>>();
        let mut kept = taken
            .iter()
            .copied()
            .filter(|contig| 2 * self.outside(*contig, offered, &rivals, &held, own) <= self.count)
            .collect::<Vec<_>>();
        if size >= 3 {
            let pool = held
                .values()
                .flatten()
                .copied()
                .filter(|contig| !primary.contains(contig))
                .collect::<BTreeSet<_>>();
            kept.extend(pool.into_iter().filter(|contig| {
                let inside = held.values().filter(|bin| bin.contains(contig)).count();
                let others = offered.holders[contig]
                    .iter()
                    .flat_map(|holder| self.voters[*holder].iter().copied())
                    .filter(|voter| !own.contains(voter) && !held.contains_key(voter))
                    .collect::<HashSet<_>>()
                    .len();
                2 * inside >= size + others
            }));
        }
        let kept = sorted(kept);
        (kept != taken).then_some(kept)
    }

    fn members(
        &self,
        judge: &Judge,
        taken: &[usize],
        primary: &HashSet<usize>,
        rivals: &HashMap<usize, Rival>,
        own: &BTreeSet<usize>,
    ) -> BTreeMap<usize, usize> {
        let size = judge.bases(taken.iter().copied());
        let mut best = BTreeMap::<usize, (usize, usize)>::new();
        for (rival, held) in rivals.iter().filter(|(_, held)| held.passes) {
            let shared = judge.bases(held.left.iter().copied().filter(|c| primary.contains(c)));
            if 2 * shared <= size.max(judge.bases(held.left.iter().copied())) {
                continue;
            }
            for voter in self.voters[*rival]
                .iter()
                .filter(|voter| !own.contains(voter))
            {
                let entry = best.entry(*voter).or_insert((*rival, shared));
                if (shared, std::cmp::Reverse(*rival)) > (entry.1, std::cmp::Reverse(entry.0)) {
                    *entry = (*rival, shared);
                }
            }
        }
        best.into_iter()
            .map(|(voter, (rival, _))| (voter, rival))
            .collect()
    }

    fn outside(
        &self,
        contig: usize,
        offered: &Offered,
        rivals: &HashMap<usize, Rival>,
        held: &BTreeMap<usize, HashSet<usize>>,
        own: &BTreeSet<usize>,
    ) -> usize {
        offered.holders[&contig]
            .iter()
            .filter(|holder| rivals.get(holder).is_some_and(|rival| rival.passes))
            .flat_map(|holder| self.voters[*holder].iter().copied())
            .filter(|voter| {
                !own.contains(voter) && held.get(voter).is_none_or(|bin| !bin.contains(&contig))
            })
            .collect::<HashSet<_>>()
            .len()
    }
}
