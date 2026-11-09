use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::quality::Scorer;
use crate::refine::owners::{credit, owners};
use crate::refine::ranking::remaining;
use crate::refine::rung::{Rung, Verdict, judge};

// Only a shard contig that fills a marker a claim lacks is mostly host. A shard that reports is a
// genome. A fused bin feeds one claim per genome and none holds most of it, so all are offered.
pub fn fold_back(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    worth: f64,
    reported: Rung,
    dissolved: &[(usize, Vec<usize>)],
    promoted: &mut [Vec<usize>],
) -> Vec<usize> {
    let claimed = promoted.iter().flatten().copied().collect::<HashSet<_>>();
    let owner = owners(promoted.iter().enumerate());
    let mut folded = Vec::new();
    for (_, contigs) in dissolved {
        let mut left = remaining(contigs, &claimed);
        if left.is_empty() || judge(features, quality, &left, reported) == Verdict::Adopt {
            continue;
        }
        let mut heirs = heirs(contigs, &owner, features, features.bin_size(&left))
            .into_iter()
            .map(|at| (at, quality.score(&promoted[at]).score(worth)))
            .collect::<Vec<_>>();
        if heirs.is_empty() {
            continue;
        }
        left.retain(|contig| !quality.features(&[*contig]).is_empty());
        left.sort_unstable_by_key(|contig| (Reverse(features.length(*contig)), *contig));
        for contig in left {
            let mut best = None;
            for (slot, (at, current)) in heirs.iter().enumerate() {
                let claim = &mut promoted[*at];
                claim.push(contig);
                let trial = quality.score(claim).score(worth);
                claim.pop();
                let gain = trial - current;
                if gain > best.map_or(0.0, |(_, best_gain, _)| best_gain) {
                    best = Some((slot, gain, trial));
                }
            }
            let Some((slot, _, trial)) = best else {
                continue;
            };
            let (at, current) = &mut heirs[slot];
            promoted[*at].push(contig);
            *current = trial;
            folded.push(contig);
        }
        for (at, _) in heirs {
            promoted[at].sort_unstable();
        }
    }
    folded
}

// A claim that took less of the bin than its shard is a stray, and the shard is another genome.
fn heirs(
    members: &[usize],
    owner: &HashMap<usize, usize>,
    features: &ContigFeatures,
    shard: usize,
) -> Vec<usize> {
    let mut held = Vec::new();
    for contig in members {
        if let Some(at) = owner.get(contig) {
            credit(&mut held, *at, features.length(*contig));
        }
    }
    held.retain(|(_, bases)| *bases > shard);
    held.sort_unstable_by_key(|(at, bases)| (Reverse(*bases), *at));
    held.into_iter().map(|(at, _)| at).collect()
}
