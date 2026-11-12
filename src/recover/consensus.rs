use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::recover::ladder::{Judge, Offered};
use crate::refine::ranking::{remaining, sorted};
use crate::tuning::FAMILY_ALPHA;

// Markers read a clean part of a genome as whole, so a winner's contigs are settled by which
// bins the other binners keep them in instead.
pub(crate) struct Ballots {
    voters: Vec<BTreeSet<usize>>,
    count: usize,
}

// The run's own bins vote first. Binners shred a host or a eukaryote alike, so their agreement on
// a piece does not outweigh the one partition built on this run's depth.
const HOME: usize = 0;

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

    pub(crate) fn several(&self) -> bool {
        self.count >= 2
    }

    pub(crate) fn support(&self, at: usize) -> usize {
        self.voters[at].len()
    }

    pub(crate) fn home(&self, at: usize) -> bool {
        self.voters[at].contains(&HOME)
    }

    pub(crate) fn agreement(&self, judge: &Judge, offered: &Offered, at: usize) -> usize {
        let size = judge.bases(offered.found[at].iter().copied());
        shared(judge, offered, &offered.found[at], |_| true)
            .into_iter()
            .filter(|(rival, held)| {
                2 * held > size.max(judge.bases(offered.found[*rival].iter().copied()))
            })
            .flat_map(|(rival, _)| self.voters[rival].iter().copied())
            .collect::<HashSet<_>>()
            .len()
    }

    pub(crate) fn resolve(
        &self,
        judge: &Judge,
        offered: &Offered,
        at: usize,
        taken: &[usize],
        claimed: &HashSet<usize>,
    ) -> Option<Vec<usize>> {
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
        let members = self.members(judge, taken, &primary, &rivals, own, true);
        if members.is_empty() {
            return None;
        }
        let size = own.len() + members.len();
        let held = members
            .iter()
            .map(|(voter, rival)| (*voter, rivals[rival].left.iter().copied().collect()))
            .collect::<BTreeMap<usize, HashSet<usize>>>();
        let voting = voting(judge, taken, &rivals);
        let mut bulk = held.clone();
        for (voter, rival) in self.members(judge, taken, &primary, &rivals, own, false) {
            bulk.entry(voter)
                .or_insert_with(|| rivals[&rival].left.iter().copied().collect());
        }
        let mut kept = taken
            .iter()
            .copied()
            .filter(|contig| {
                let inside = own.len() + held.values().filter(|bin| bin.contains(contig)).count();
                inside >= self.outside(*contig, offered, &voting, &rivals, &bulk, own)
            })
            .collect::<Vec<_>>();
        // A vote that keeps none of a passing winner claims nothing, so the winner would pop again.
        if kept.is_empty() {
            return None;
        }
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
        let kept = match std::env::var_os("ROSELLA_VOTE_TIER") {
            Some(_) => held_tier(judge, taken, sorted(kept)),
            None => sorted(kept),
        };
        (kept != taken).then_some(kept)
    }

    // Markers and depth cannot tell a fused bin under the floor from a whole genome, so the cut
    // goes by which binners offer a small piece of it apart. They shred a eukaryote alike.
    pub(crate) fn pieces(
        &self,
        judge: &Judge,
        offered: &Offered,
        taken: &[usize],
        claimed: &HashSet<usize>,
        free: &[usize],
        open: impl Fn(usize) -> bool,
    ) -> Vec<(usize, Vec<usize>)> {
        let size = judge.bases(taken.iter().copied());
        let primary = taken.iter().copied().collect::<HashSet<_>>();
        let inside = shared(judge, offered, taken, |_| true);
        let mut cuts = inside
            .iter()
            .filter(|(rival, _)| open(**rival))
            .filter(|(rival, held)| 2 * **held > free[**rival] && 2 * free[**rival] < size)
            .filter_map(|(rival, held)| {
                let left = remaining(&offered.found[*rival], claimed);
                let within = shared(judge, offered, &left, |_| true);
                let within_winner = shared(judge, offered, &left, |at| primary.contains(&at));
                let rest = size - held;
                let (mut joins, mut splits) = (BTreeSet::new(), BTreeSet::new());
                for (holder, part) in &within {
                    let outside = inside.get(holder).copied().unwrap_or(0)
                        - within_winner.get(holder).copied().unwrap_or(0);
                    if 2 * part > free[*rival] && 2 * outside > rest {
                        joins.extend(self.voters[*holder].iter().copied());
                    } else if 2 * part > free[*rival].max(free[*holder]) {
                        splits.extend(self.voters[*holder].iter().copied());
                    }
                }
                let splits = splits.difference(&joins).count();
                (splits > joins.len() && !joins.contains(&HOME))
                    .then(|| (splits - joins.len(), free[*rival], *rival, left))
            })
            .collect::<Vec<_>>();
        cuts.sort_by_key(|(margin, bases, rival, _)| {
            (
                std::cmp::Reverse(*margin),
                std::cmp::Reverse(*bases),
                *rival,
            )
        });
        let mut used = HashSet::new();
        cuts.into_iter()
            .filter(|(.., left)| {
                let apart = left.iter().all(|contig| !used.contains(contig));
                if apart {
                    used.extend(left.iter().copied());
                }
                apart
            })
            .map(|(_, _, rival, left)| (rival, left))
            .collect()
    }

    fn members(
        &self,
        judge: &Judge,
        taken: &[usize],
        primary: &HashSet<usize>,
        rivals: &HashMap<usize, Rival>,
        own: &BTreeSet<usize>,
        passing: bool,
    ) -> BTreeMap<usize, usize> {
        let size = judge.bases(taken.iter().copied());
        let mut best = BTreeMap::<usize, (usize, usize)>::new();
        for (rival, held) in rivals.iter().filter(|(_, held)| held.passes == passing) {
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
        voting: &HashSet<usize>,
        rivals: &HashMap<usize, Rival>,
        bulk: &BTreeMap<usize, HashSet<usize>>,
        own: &BTreeSet<usize>,
    ) -> usize {
        // A binner that shreds the winner cut it along depth, so its pieces sit apart and say
        // nothing. A failing bin votes only beside its binner's one bin holding most of the winner.
        offered.holders[&contig]
            .iter()
            .filter(|holder| voting.contains(holder))
            .flat_map(|holder| {
                let passes = rivals[holder].passes;
                self.voters[*holder]
                    .iter()
                    .copied()
                    .filter(move |voter| passes || bulk.contains_key(voter))
            })
            .filter(|voter| {
                !own.contains(voter) && bulk.get(voter).is_none_or(|bin| !bin.contains(&contig))
            })
            .collect::<HashSet<_>>()
            .len()
    }
}

// Markers cannot see marker-free contigs, so binners vote on those. Marker contigs are the
// markers' own evidence, so the vote never takes a winner out of the counted tier.
fn held_tier(judge: &Judge, taken: &[usize], kept: Vec<usize>) -> Vec<usize> {
    let accept = judge.bars.at(0).completeness;
    let counted = |contigs: &[usize]| {
        let quality = judge.score(contigs);
        quality.completeness >= accept && quality.contamination <= judge.bars.tier()
    };
    if !counted(taken) || counted(&kept) {
        return kept;
    }
    let inside = kept.iter().copied().collect::<HashSet<_>>();
    let carriers = taken
        .iter()
        .copied()
        .filter(|at| !inside.contains(at) && judge.markers.carries(judge.contigs[*at]));
    let back = sorted(kept.iter().copied().chain(carriers));
    match counted(&back) {
        true => back,
        false => taken.to_vec(),
    }
}

fn shared(
    judge: &Judge,
    offered: &Offered,
    contigs: &[usize],
    counted: impl Fn(usize) -> bool,
) -> BTreeMap<usize, usize> {
    let mut held = BTreeMap::new();
    for contig in contigs.iter().copied().filter(|contig| counted(*contig)) {
        let bases = judge.bases([contig]);
        for holder in &offered.holders[&contig] {
            *held.entry(*holder).or_default() += bases;
        }
    }
    held
}

// A marker-free element's own bin can never pass, so a failing bin votes too, but only when its
// piece of the winner sits apart in depth: a binner's offcut of the same genome reads level.
fn voting(judge: &Judge, taken: &[usize], rivals: &HashMap<usize, Rival>) -> HashSet<usize> {
    let failing = rivals.values().filter(|rival| !rival.passes).count();
    let profile = (failing > 0).then(|| Profile::new(judge, taken));
    rivals
        .par_iter()
        .filter(|(id, rival)| {
            rival.passes
                || profile
                    .as_ref()
                    .is_some_and(|profile| profile.apart(taken, &rival.left, failing, **id as u64))
        })
        .map(|(id, _)| *id)
        .collect()
}

struct Profile {
    rows: Vec<Vec<f64>>,
    total: Vec<f64>,
}

impl Profile {
    fn new(judge: &Judge, taken: &[usize]) -> Self {
        let features = judge.depth;
        let rows = taken
            .iter()
            .map(|at| {
                features
                    .coverage_row(judge.contigs[*at])
                    .as_chunks::<2>()
                    .0
                    .iter()
                    .map(|sample| sample[0].ln_1p())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut total = vec![0.0; features.n_samples()];
        for row in &rows {
            for (slot, value) in total.iter_mut().zip(row) {
                *slot += value;
            }
        }
        Self { rows, total }
    }

    fn gap(&self, picked: &[usize]) -> f64 {
        let mut sum = vec![0.0; self.total.len()];
        for row in picked {
            for (slot, value) in sum.iter_mut().zip(&self.rows[*row]) {
                *slot += value;
            }
        }
        let inside = picked.len() as f64;
        let rest = (self.rows.len() - picked.len()) as f64;
        sum.iter()
            .zip(&self.total)
            .map(|(part, whole)| (part / inside - (whole - part) / rest).powi(2))
            .sum()
    }

    // Same-sized random draws from the winner are the null, sized so that one draw resolves the
    // corrected level, as the splitter's dip test is.
    fn apart(&self, taken: &[usize], left: &[usize], tests: usize, seed: u64) -> bool {
        let left = left.iter().collect::<HashSet<_>>();
        let picked = taken
            .iter()
            .enumerate()
            .filter(|(_, at)| left.contains(at))
            .map(|(row, _)| row)
            .collect::<Vec<_>>();
        let n = self.rows.len();
        if picked.is_empty() || picked.len() == n {
            return false;
        }
        let observed = self.gap(&picked);
        let draws = (tests as f64 / FAMILY_ALPHA).ceil() as usize;
        let mut rng = StdRng::seed_from_u64(seed);
        let mut order = (0..n).collect::<Vec<_>>();
        (0..draws).all(|_| {
            for slot in 0..picked.len() {
                order.swap(slot, rng.random_range(slot..n));
            }
            self.gap(&order[..picked.len()]) < observed
        })
    }
}
