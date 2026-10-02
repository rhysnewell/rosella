use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap, HashSet};

use log::debug;
use rayon::prelude::*;

use crate::clustering::clusterer::Partitioning;
use crate::clustering::graph_partition::Partition;
use crate::quality::{Quality, Scorer};
use crate::recover::combine_report::CombineReport;
use crate::refine::ranking::{Ranked, remaining, sorted};
use crate::refine::rung::Bars;

// Codelength picks a rung about half as fine as the truth, so each arm contributes the rung
// its markers choose rather than the one codelength ranks first.
pub fn best_per_arm(ladder: Vec<Partitioning>, judge: &Judge) -> Vec<Partitioning> {
    let mut arms: Vec<((Partition, u64), Vec<Partitioning>)> = Vec::new();
    for held in ladder {
        let key = (held.arm, held.seed);
        match arms.iter_mut().find(|(seen, _)| *seen == key) {
            Some((_, entries)) => entries.push(held),
            None => arms.push((key, vec![held])),
        }
    }
    arms.into_iter()
        .map(|(_, entries)| pick_rung(entries, judge))
        .collect()
}

pub struct Judge<'a> {
    pub quality: &'a dyn Scorer,
    pub contigs: &'a [usize],
    pub bars: Bars,
}

impl Judge<'_> {
    pub(crate) fn score(&self, positions: &[usize]) -> Quality {
        let mapped = positions
            .iter()
            .map(|at| self.contigs[*at])
            .collect::<Vec<_>>();
        self.quality.score(&mapped)
    }

    fn worth(&self, positions: &[usize]) -> f64 {
        self.score(positions).score(self.bars.worth)
    }

    pub(crate) fn mapped(&self, positions: &[usize]) -> Vec<usize> {
        positions.iter().map(|at| self.contigs[*at]).collect()
    }
}

// Codelength ranks every rung about half as fine as the truth, so where an annotation exists
// the markers judge the ladder instead.
pub fn pick_rung(ladder: Vec<Partitioning>, judge: &Judge) -> Partitioning {
    let bars = crate::quality::Bars {
        completeness: judge.bars.completeness,
        contamination: judge.bars.tier(),
    };
    let mut scored = ladder
        .into_par_iter()
        .map(|held| {
            let count = held
                .cluster_map
                .par_iter()
                .filter(|(_, members)| judge.score(members).clears(bars))
                .count();
            debug!(
                "rung {} communities, {count} over the bar",
                held.cluster_map.len()
            );
            (count, held)
        })
        .collect::<Vec<_>>();
    scored.sort_by_key(|(count, _)| Reverse(*count));
    let (count, chosen) = scored.swap_remove(0);
    debug!(
        "Markers chose a {} community rung, {count} bins over the bar.",
        chosen.cluster_map.len()
    );
    chosen
}

pub(crate) fn candidates(ladder: &[Partitioning]) -> Vec<Vec<usize>> {
    let mut found = ladder
        .iter()
        .flat_map(|held| held.cluster_map.values().cloned())
        .collect::<Vec<_>>();
    found.sort_unstable();
    found.dedup();
    found
}

// Choosing one rung whole is worth almost nothing against choosing the best of both arms, so the
// bins are arbitrated one at a time instead and a rung contributes only the ones that win.
pub fn combine(
    ladder: &[Partitioning],
    judge: &Judge,
    report: Option<&CombineReport<'_>>,
) -> Partitioning {
    let _timer = crate::timing::scope("combine");
    let every = ladder
        .iter()
        .flat_map(|held| held.cluster_map.values().flatten().copied())
        .chain(ladder.iter().flat_map(|held| held.outliers.iter().copied()))
        .collect::<HashSet<_>>();
    let mut held = BinaryHeap::from(
        candidates(ladder)
            .into_par_iter()
            .map(|contigs| Ranked {
                worth: judge.worth(&contigs),
                contigs,
                extra: (),
            })
            .collect::<Vec<_>>(),
    );

    let mut cluster_map = BTreeMap::new();
    let mut claimed = HashSet::new();
    let mut seen = 0;
    while let Some(entry) = held.pop() {
        let left = remaining(&entry.contigs, &claimed);
        if let Some(report) = report {
            let verdict = if left.is_empty() {
                "spent"
            } else if left.len() < entry.contigs.len() {
                "trimmed"
            } else {
                "taken"
            };
            report.row(seen, entry.worth, verdict, &judge.mapped(&entry.contigs));
            seen += 1;
        }
        if left.is_empty() {
            continue;
        }
        if left.len() < entry.contigs.len() {
            held.push(Ranked {
                worth: judge.worth(&left),
                contigs: left,
                extra: (),
            });
            continue;
        }
        claimed.extend(left.iter().copied());
        cluster_map.insert(cluster_map.len(), left);
    }

    let outliers = sorted(every.difference(&claimed).copied());
    debug!(
        "Markers combined {} rungs into {} bins, {} contigs unclaimed.",
        ladder.len(),
        cluster_map.len(),
        outliers.len()
    );
    Partitioning {
        cluster_map,
        outliers,
        score: None,
        arm: Partition::Both,
        seed: 0,
    }
}
