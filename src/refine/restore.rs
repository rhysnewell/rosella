use std::collections::{HashMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::quality::{Scorer, Worth};
use crate::refine::rung::{Rung, Verdict, judge};

pub struct Restored {
    pub promoted: Vec<Vec<usize>>,
    pub released: Vec<usize>,
    pub bins: usize,
}

struct Thread {
    bins: Vec<usize>,
    pieces: Vec<usize>,
}

fn components(dissolved: &[(usize, Vec<usize>)], draws: &[HashSet<usize>]) -> Vec<Thread> {
    let mut holders: HashMap<usize, Vec<usize>> = HashMap::new();
    for (at, taken) in draws.iter().enumerate() {
        for bin in taken {
            holders.entry(*bin).or_default().push(at);
        }
    }
    let mut seen = HashSet::new();
    let mut found = Vec::new();
    for (bin, _) in dissolved {
        if !holders.contains_key(bin) || !seen.insert(*bin) {
            continue;
        }
        let (mut bins, mut pieces, mut queue) = (vec![*bin], HashSet::new(), vec![*bin]);
        while let Some(bin) = queue.pop() {
            for at in holders.get(&bin).into_iter().flatten() {
                if !pieces.insert(*at) {
                    continue;
                }
                for reached in &draws[*at] {
                    if seen.insert(*reached) {
                        bins.push(*reached);
                        queue.push(*reached);
                    }
                }
            }
        }
        let mut pieces = pieces.into_iter().collect::<Vec<_>>();
        pieces.sort_unstable();
        bins.sort_unstable();
        found.push(Thread { bins, pieces });
    }
    found
}

/// Finished genomes first, then everything the run would still report, then the best single bin.
/// The reporting tier is what makes breaking up clusters nothing counts free, and it sits under
/// the accept bar so two halves of a genome never outvote the genome.
fn state(
    held: &Judge,
    worth: Worth,
    bins: &[Vec<usize>],
) -> (usize, usize, f64) {
    let over = |bar: Rung| {
        bins.iter()
            .filter(|contigs| {
                !contigs.is_empty()
                    && judge(held.features, held.quality, contigs, bar) == Verdict::Adopt
            })
            .count()
    };
    let best = bins
        .iter()
        .filter(|contigs| !contigs.is_empty())
        .map(|contigs| held.quality.score(contigs).score(worth))
        .fold(f64::NEG_INFINITY, f64::max);
    (over(held.accept), over(held.countable), best)
}

pub struct Judge<'a> {
    pub features: &'a ContigFeatures<'a>,
    pub quality: &'a dyn Scorer,
    pub countable: Rung,
    pub accept: Rung,
}

/// The pool hands a dissolved bin its leftovers, never itself. A piece can draw from several
/// bins, so the unit that can be reverted is the whole connected run of bins and pieces, and it
/// is kept only when the pool's arrangement of it beats the bins it was made from.
pub fn restore(
    held: &Judge,
    worth: Worth,
    dissolved: &[(usize, Vec<usize>)],
    promoted: Vec<Vec<usize>>,
) -> Restored {
    let origin = dissolved
        .iter()
        .flat_map(|(bin, contigs)| contigs.iter().map(|contig| (*contig, *bin)))
        .collect::<HashMap<_, _>>();
    let bin_of = dissolved
        .iter()
        .map(|(bin, contigs)| (*bin, contigs))
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

    let mut dropped = HashSet::new();
    let mut bins = 0;
    for thread in components(dissolved, &draws) {
        let claimed = thread
            .pieces
            .iter()
            .flat_map(|at| promoted[*at].iter().copied())
            .collect::<HashSet<_>>();
        let mut keep = thread
            .pieces
            .iter()
            .map(|at| promoted[*at].clone())
            .collect::<Vec<_>>();
        let mut revert = Vec::new();
        for bin in &thread.bins {
            let Some(contigs) = bin_of.get(bin) else {
                continue;
            };
            keep.push(
                contigs
                    .iter()
                    .copied()
                    .filter(|contig| !claimed.contains(contig))
                    .collect(),
            );
            revert.push((*contigs).clone());
        }
        if state(held, worth, &revert) <= state(held, worth, &keep) {
            continue;
        }
        dropped.extend(thread.pieces);
        bins += revert.len();
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
