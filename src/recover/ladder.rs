use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};

use log::debug;

use crate::clustering::clusterer::Partitioning;
use crate::clustering::graph_partition::Partition;
use crate::quality::{Quality, Scorer};
use crate::refine::rung::Bars;
use crate::refine::select::{Ranked, remaining, sorted};

/// Codelength picks a rung about half as fine as the truth, so each arm contributes the rung
/// its markers choose rather than the one codelength ranks first.
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
    fn score(&self, positions: &[usize]) -> Quality {
        let mapped = positions
            .iter()
            .map(|at| self.contigs[*at])
            .collect::<Vec<_>>();
        self.quality.score(&mapped)
    }

    fn worth(&self, positions: &[usize]) -> f64 {
        self.score(positions).score(self.bars.worth)
    }
}

/// Codelength ranks every rung about half as fine as the truth, so where an annotation exists
/// the markers judge the ladder instead.
pub fn pick_rung(ladder: Vec<Partitioning>, judge: &Judge) -> Partitioning {
    let bar = judge.bars.completeness;
    let tier = judge.bars.tier();
    let passing = |held: &Partitioning| {
        held.cluster_map
            .values()
            .map(|members| judge.score(&sorted(members.iter().copied())))
            .filter(|held| held.completeness >= bar && held.contamination <= tier)
            .count()
    };
    let mut scored = ladder
        .into_iter()
        .map(|held| {
            let count = passing(&held);
            debug!(
                "rung {} communities, {count} over the bar",
                held.cluster_map.len()
            );
            (count as f64, held)
        })
        .collect::<Vec<_>>();
    scored.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(Ordering::Equal));
    let (value, chosen) = scored.swap_remove(0);
    debug!(
        "Markers chose a {} community rung, {value:.0} bins over the bar.",
        chosen.cluster_map.len()
    );
    chosen
}

fn candidates(ladder: &[Partitioning]) -> Vec<Vec<usize>> {
    let mut found = ladder
        .iter()
        .flat_map(|held| {
            held.cluster_map
                .values()
                .map(|members| sorted(members.iter().copied()))
        })
        .collect::<Vec<_>>();
    found.sort_unstable();
    found.dedup();
    found
}

/// Choosing one rung whole is worth almost nothing against choosing the best of both arms, so the
/// bins are arbitrated one at a time instead and a rung contributes only the ones that win.
pub fn combine(ladder: Vec<Partitioning>, judge: &Judge) -> Partitioning {
    let _timer = crate::timing::scope("combine");
    let every = ladder
        .iter()
        .flat_map(|held| held.cluster_map.values().flatten().copied())
        .chain(ladder.iter().flat_map(|held| held.outliers.iter().copied()))
        .collect::<HashSet<_>>();
    let mut held = candidates(&ladder)
        .into_iter()
        .map(|contigs| Ranked {
            worth: judge.worth(&contigs),
            contigs,
            extra: (),
        })
        .collect::<BinaryHeap<_>>();

    let mut cluster_map: HashMap<usize, HashSet<usize>> = HashMap::new();
    let mut claimed = HashSet::new();
    while let Some(entry) = held.pop() {
        let left = remaining(&entry.contigs, &claimed);
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
        cluster_map.insert(cluster_map.len(), left.into_iter().collect());
    }

    let outliers = every.difference(&claimed).copied().collect::<HashSet<_>>();
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
