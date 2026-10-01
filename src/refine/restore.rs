use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::quality::Scorer;
use crate::refine::owners::owners;
use crate::refine::rung::{Rung, Verdict, verdict};

// A genome with duplicated marker families repeats sequence at this rate; below the floor the
// bin is more often a chimera of unrelated genomes, which repeats nothing either.
const PURE_BAR: f64 = 0.05;
const PURE_FLOOR: f64 = 0.04;

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
    // A genome's own duplicated single copy families read as contamination, so k-mers decide
    // where contamination is the only objection. An incomplete bin is the pool's job.
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

    // The best single bin decides, and the counts only break its ties. Counting first rewards
    // cutting a genome in two, since both halves report, where worth never does.
    fn state(&self, worth: f64, bins: &[Vec<usize>]) -> (f64, usize, usize) {
        let (mut best, mut accepted, mut reported) = (f64::NEG_INFINITY, 0, 0);
        for contigs in bins.iter().filter(|contigs| !contigs.is_empty()) {
            let quality = self.quality.score(contigs);
            let bases = self.features.bin_size(contigs);
            let adopted = |rung: Rung| usize::from(verdict(bases, quality, rung) == Verdict::Adopt);
            best = best.max(quality.score(worth));
            accepted += adopted(self.accept);
            reported += adopted(self.reported);
        }
        (best, accepted, reported)
    }
}

fn better(left: (f64, usize, usize), right: (f64, usize, usize)) -> bool {
    left.0
        .total_cmp(&right.0)
        .then(left.1.cmp(&right.1))
        .then(left.2.cmp(&right.2))
        == Ordering::Greater
}

// Dropping a pool piece hurts no other bin, since its contigs go back where they came from. So one
// dissolved bin is weighed against the pieces holding its contigs, scoring the bins they touch.
pub fn restore(
    held: &Judge,
    worth: f64,
    dissolved: &[(usize, Vec<usize>)],
    promoted: Vec<Vec<usize>>,
) -> Restored {
    let origin = owners(dissolved.iter().map(|(bin, contigs)| (*bin, contigs)));
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

    let mut order = dissolved
        .iter()
        .map(|entry| (held.worth_of(worth, &entry.1), entry))
        .collect::<Vec<_>>();
    order.sort_by(|(ours, (left, _)), (theirs, (right, _))| {
        theirs.total_cmp(ours).then(left.cmp(right))
    });

    // Rebuilding the live set for every dissolved bin was quadratic in the promoted contigs.
    let mut live = HashMap::<usize, usize>::new();
    for contig in promoted.iter().flatten() {
        *live.entry(*contig).or_default() += 1;
    }
    let mut dropped = HashSet::new();
    let mut bins = 0;
    for (_, (bin, contigs)) in order {
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
            release(&pieces, &promoted, &mut live);
            dropped.extend(pieces);
            bins += 1;
            continue;
        }

        let touched = pieces
            .iter()
            .flat_map(|at| draws[*at].iter().copied())
            .collect::<HashSet<_>>();

        let mut freed = HashMap::<usize, usize>::new();
        for contig in pieces.iter().flat_map(|at| &promoted[*at]) {
            *freed.entry(*contig).or_default() += 1;
        }
        let held_now = |contig: &usize| live.get(contig).is_some_and(|count| *count > 0);
        let held_after = |contig: &usize| {
            live.get(contig).copied().unwrap_or(0) > freed.get(contig).copied().unwrap_or(0)
        };

        let mut keep = pieces
            .iter()
            .map(|at| promoted[*at].clone())
            .collect::<Vec<_>>();
        let mut revert = vec![contigs.clone()];
        for (other, theirs) in dissolved
            .iter()
            .filter(|(other, _)| touched.contains(other))
        {
            keep.push(unheld(theirs, held_now));
            if other != bin {
                revert.push(unheld(theirs, held_after));
            }
        }
        if !better(held.state(worth, &revert), held.state(worth, &keep)) {
            continue;
        }
        release(&pieces, &promoted, &mut live);
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

fn release(pieces: &[usize], promoted: &[Vec<usize>], live: &mut HashMap<usize, usize>) {
    for contig in pieces.iter().flat_map(|at| &promoted[*at]) {
        if let Some(count) = live.get_mut(contig) {
            *count -= 1;
        }
    }
}

fn unheld(contigs: &[usize], held: impl Fn(&usize) -> bool) -> Vec<usize> {
    contigs
        .iter()
        .copied()
        .filter(|contig| !held(contig))
        .collect()
}
