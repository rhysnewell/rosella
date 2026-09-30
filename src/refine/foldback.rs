use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::quality::Scorer;
use crate::refine::rung::{Rung, Verdict, judge};
use crate::refine::select::remaining;

// The shard a claim leaves of its bin holds more foreign bases than host, and only a contig that
// fills a marker the claim lacks is mostly host. A leftover that reports as a bin is a genome.
pub fn fold_back(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    worth: f64,
    reported: Rung,
    dissolved: &[(usize, Vec<usize>)],
    promoted: &mut [Vec<usize>],
) -> Vec<usize> {
    let claimed = promoted.iter().flatten().copied().collect::<HashSet<_>>();
    let owner = promoted
        .iter()
        .enumerate()
        .flat_map(|(at, contigs)| contigs.iter().map(move |contig| (*contig, at)))
        .collect::<HashMap<_, _>>();
    let mut folded = Vec::new();
    for (_, contigs) in dissolved {
        let mut left = remaining(contigs, &claimed);
        if left.is_empty() || judge(features, quality, &left, reported) == Verdict::Adopt {
            continue;
        }
        let mut taken = HashMap::<usize, usize>::new();
        for contig in contigs {
            if let Some(at) = owner.get(contig) {
                *taken.entry(*at).or_default() += features.length(*contig);
            }
        }
        let Some((at, held)) = taken
            .into_iter()
            .max_by_key(|(at, held)| (*held, Reverse(*at)))
        else {
            continue;
        };
        if 2 * held <= features.bin_size(contigs) {
            continue;
        }
        left.sort_unstable_by_key(|contig| (Reverse(features.length(*contig)), *contig));
        let claim = &mut promoted[at];
        let mut current = quality.score(claim).score(worth);
        for contig in left {
            claim.push(contig);
            let trial = quality.score(claim).score(worth);
            if trial <= current {
                claim.pop();
                continue;
            }
            current = trial;
            folded.push(contig);
        }
        claim.sort_unstable();
    }
    folded
}
