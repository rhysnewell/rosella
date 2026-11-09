use std::cmp::Reverse;
use std::collections::HashSet;

use crate::embedding::features::ContigFeatures;
use crate::quality::Scorer;
use crate::refine::owners::{heir, owners};
use crate::refine::ranking::remaining;
use crate::refine::rung::{Rung, Verdict, judge};

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
    let owner = owners(promoted.iter().enumerate());
    let mut folded = Vec::new();
    for (_, contigs) in dissolved {
        let mut left = remaining(contigs, &claimed);
        if left.is_empty() || judge(features, quality, &left, reported) == Verdict::Adopt {
            continue;
        }
        let Some((at, held)) = heir(contigs, &owner, |contig| features.length(contig)) else {
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
