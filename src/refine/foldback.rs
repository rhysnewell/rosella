use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::embedding::metrics::{AggregateMetric, Point};
use crate::quality::Scorer;
use crate::refine::bin_stats::{Weighted, centroid};
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
    let metric = AggregateMetric::new(
        features.n_samples() * 2,
        features.distance_settings().clone(),
    );
    let mut folded = Vec::new();
    for (_, contigs) in dissolved {
        let mut left = remaining(contigs, &claimed);
        if left.is_empty() || judge(features, quality, &left, reported) == Verdict::Adopt {
            continue;
        }
        let shard = features.bin_size(&left);
        let mut heirs = heirs(contigs, &owner, features)
            .into_iter()
            .map(|(at, took)| Heir {
                at,
                current: quality.score(&promoted[at]).score(worth),
                centre: (took <= shard).then(|| centroid(features, &promoted[at])),
            })
            .collect::<Vec<_>>();
        if heirs.is_empty() {
            continue;
        }
        let rest = Weighted::of(features, &left);
        left.retain(|contig| !quality.features(&[*contig]).is_empty());
        left.sort_unstable_by_key(|contig| (Reverse(features.length(*contig)), *contig));
        for contig in left {
            let point = features.point(contig);
            let home = rest
                .without(features, contig)
                .map_or(f64::INFINITY, |rest| metric.distance(&point, &rest));
            let mut best = None;
            for (slot, heir) in heirs.iter().enumerate() {
                if heir
                    .centre
                    .as_ref()
                    .is_some_and(|centre| metric.distance(&point, centre) > home)
                {
                    continue;
                }
                let claim = &mut promoted[heir.at];
                claim.push(contig);
                let trial = quality.score(claim).score(worth);
                claim.pop();
                let gain = trial - heir.current;
                if gain > best.map_or(0.0, |(_, best_gain, _)| best_gain) {
                    best = Some((slot, gain, trial));
                }
            }
            let Some((slot, _, trial)) = best else {
                continue;
            };
            let heir = &mut heirs[slot];
            promoted[heir.at].push(contig);
            heir.current = trial;
            folded.push(contig);
        }
        for heir in heirs {
            promoted[heir.at].sort_unstable();
        }
    }
    folded
}

// A claim that took less of the bin than its shard holds may face another genome's homologue, so
// it must sit nearer the contig than the shard does. A larger claim's shard is its own remnant.
struct Heir {
    at: usize,
    current: f64,
    centre: Option<Point>,
}

// The claim holding most of the bin is offered first, so it wins a tie.
fn heirs(
    members: &[usize],
    owner: &HashMap<usize, usize>,
    features: &ContigFeatures,
) -> Vec<(usize, usize)> {
    let mut held = Vec::new();
    for contig in members {
        if let Some(at) = owner.get(contig) {
            credit(&mut held, *at, features.length(*contig));
        }
    }
    held.sort_unstable_by_key(|(at, bases)| (Reverse(*bases), *at));
    held
}
