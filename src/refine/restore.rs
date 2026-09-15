use std::collections::{HashMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::quality::Scorer;
use crate::refine::rung::{Rung, Verdict, judge};

fn clears(features: &ContigFeatures, quality: &dyn Scorer, contigs: &[usize], bar: Rung) -> bool {
    !contigs.is_empty() && judge(features, quality, contigs, bar) == Verdict::Adopt
}

pub struct Restored {
    pub promoted: Vec<Vec<usize>>,
    pub released: Vec<usize>,
    pub bins: usize,
}

/// The pool hands a dissolved bin its leftovers, never itself, so a bin it takes apart into
/// pieces that all miss the bar is a genome lost to no one. Put those back whole.
pub fn restore(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    dissolved: &[(usize, Vec<usize>)],
    promoted: Vec<Vec<usize>>,
    bar: Rung,
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

    let mut dropped = HashSet::new();
    let mut bins = 0;
    for (bin, contigs) in dissolved {
        if !clears(features, quality, contigs, bar) {
            continue;
        }
        let pieces = draws
            .iter()
            .enumerate()
            .filter(|(_, taken)| taken.contains(bin))
            .map(|(at, _)| at)
            .collect::<Vec<_>>();
        // A piece drawing from two bins cannot be unpicked for one of them without orphaning
        // the other's contigs, and a piece already promised to another bin is spoken for.
        if pieces
            .iter()
            .any(|at| draws[*at].len() > 1 || dropped.contains(at))
        {
            continue;
        }
        let claimed = pieces
            .iter()
            .flat_map(|at| promoted[*at].iter().copied())
            .collect::<HashSet<_>>();
        let remnant = contigs
            .iter()
            .copied()
            .filter(|contig| !claimed.contains(contig))
            .collect::<Vec<_>>();
        let kept = pieces
            .iter()
            .filter(|at| clears(features, quality, &promoted[**at], bar))
            .count()
            + usize::from(clears(features, quality, &remnant, bar));
        if kept > 0 {
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
