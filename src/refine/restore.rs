use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::quality::Scorer;
use crate::refine::rung::{Rung, Verdict, judge};
use crate::refine::select::remaining;

/// A genome with duplicated marker families repeats sequence at this rate; below the floor the
/// bin is more often a chimera of unrelated genomes, which repeats nothing either.
pub const PURE_BAR: f64 = 0.05;
pub const PURE_FLOOR: f64 = 0.04;

pub struct Restored {
    pub promoted: Vec<Vec<usize>>,
    pub released: Vec<usize>,
    pub bins: usize,
}

pub struct Judge<'a> {
    pub features: &'a ContigFeatures<'a>,
    pub quality: &'a dyn Scorer,
    pub reported: Rung,
    pub accept: Rung,
}

impl Judge<'_> {
    /// A genome's own duplicated single copy families read as contamination, so k-mers decide
    /// where contamination is the only objection. An incomplete bin is the pool's job.
    fn one_organism(&self, contigs: &[usize]) -> bool {
        let quality = self.quality.score(contigs);
        if quality.completeness < self.accept.completeness
            || quality.contamination <= self.accept.contamination
        {
            return false;
        }
        self.features
            .sketches()
            .and_then(|sketches| sketches.duplication(contigs))
            .is_some_and(|share| (PURE_FLOOR..=PURE_BAR).contains(&share))
    }

    fn worth_of(&self, worth: f64, contigs: &[usize]) -> f64 {
        self.quality.score(contigs).score(worth)
    }

    /// The best single bin decides, and the counts only break its ties. Counting first rewards
    /// cutting a genome in two, since both halves report, where worth never does.
    fn state(&self, worth: f64, bins: &[Vec<usize>]) -> (f64, usize, usize) {
        let over = |bar: Rung| {
            bins.iter()
                .filter(|contigs| {
                    !contigs.is_empty()
                        && judge(self.features, self.quality, contigs, bar) == Verdict::Adopt
                })
                .count()
        };
        let best = bins
            .iter()
            .filter(|contigs| !contigs.is_empty())
            .map(|contigs| self.worth_of(worth, contigs))
            .fold(f64::NEG_INFINITY, f64::max);
        (best, over(self.accept), over(self.reported))
    }
}

fn better(left: (f64, usize, usize), right: (f64, usize, usize)) -> bool {
    left.0
        .total_cmp(&right.0)
        .then(left.1.cmp(&right.1))
        .then(left.2.cmp(&right.2))
        == Ordering::Greater
}

/// The pool hands a dissolved bin its leftovers, never itself. Dropping a piece hurts no other
/// bin, since every contig in it goes back where it came from, so the unit weighed here is one
/// bin against the pieces holding its contigs, with the other bins those pieces touch scored
/// either way round.
pub fn restore(
    held: &Judge,
    worth: f64,
    dissolved: &[(usize, Vec<usize>)],
    promoted: Vec<Vec<usize>>,
) -> Restored {
    let origin = dissolved
        .iter()
        .flat_map(|(bin, contigs)| contigs.iter().map(|contig| (*contig, *bin)))
        .collect::<HashMap<_, _>>();
    let draws = promoted
        .iter()
        .map(|contigs| {
            contigs
                .iter()
                .filter_map(|contig| origin.get(contig))
                .copied()
                .collect::<HashSet<_>>()
        })
        .collect::<Vec<_>>();

    let mut order = dissolved.iter().collect::<Vec<_>>();
    order.sort_by(|(left, ours), (right, theirs)| {
        held.worth_of(worth, theirs)
            .total_cmp(&held.worth_of(worth, ours))
            .then(left.cmp(right))
    });

    let mut dropped = HashSet::new();
    let mut bins = 0;
    for (bin, contigs) in order {
        let pieces = draws
            .iter()
            .enumerate()
            .filter(|(at, taken)| taken.contains(bin) && !dropped.contains(at))
            .map(|(at, _)| at)
            .collect::<Vec<_>>();
        if pieces.is_empty() {
            continue;
        }
        if held.one_organism(contigs) {
            dropped.extend(pieces);
            bins += 1;
            continue;
        }

        let touched = pieces
            .iter()
            .flat_map(|at| draws[*at].iter().copied())
            .collect::<HashSet<_>>();

        let live = |without: &HashSet<usize>| {
            draws
                .iter()
                .enumerate()
                .filter(|(at, _)| !dropped.contains(at) && !without.contains(at))
                .flat_map(|(at, _)| promoted[at].iter().copied())
                .collect::<HashSet<_>>()
        };
        let now = live(&HashSet::new());
        let after = live(&pieces.iter().copied().collect());

        let mut keep = pieces
            .iter()
            .map(|at| promoted[*at].clone())
            .collect::<Vec<_>>();
        let mut revert = vec![contigs.clone()];
        for (other, theirs) in dissolved.iter().filter(|(other, _)| touched.contains(other)) {
            keep.push(remaining(theirs, &now));
            if other != bin {
                revert.push(remaining(theirs, &after));
            }
        }
        if !better(held.state(worth, &revert), held.state(worth, &keep)) {
            continue;
        }
        dropped.extend(pieces);
        bins += 1;
    }

    let mut released = Vec::new();
    let promoted = promoted
        .into_iter()
        .enumerate()
        .filter_map(|(at, contigs)| match dropped.contains(&at) {
            true => {
                released.extend(contigs);
                None
            }
            false => Some(contigs),
        })
        .collect();
    Restored {
        promoted,
        released,
        bins,
    }
}
