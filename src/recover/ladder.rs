use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};

use log::{debug, info};

use crate::clustering::clusterer::Partitioning;
use crate::clustering::graph_partition::Partition;
use crate::quality::{Quality, Scorer, Worth};
use crate::refine::select::remaining;

/// The tier a recovered genome is counted at, not the accept bar, because the rung is being
/// judged on how many genomes it would yield rather than on what the pool will take.
const TIER_CONTAMINATION: f64 = 10.0;

/// The graph objective picks a rung about half as fine as the truth, so each arm contributes the
/// rung its markers choose rather than the one the objective ranks first.
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
    pub worth: Worth,
    pub completeness: f64,
    pub contamination: f64,
    pub bar: bool,
    pub size_tie: bool,
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
        self.score(positions).score(self.worth)
    }

    fn admits(&self, positions: &[usize]) -> bool {
        if !self.bar {
            return true;
        }
        let held = self.score(positions);
        held.completeness >= self.quality.completeness_bar(self.completeness)
            && held.contamination <= self.contamination
    }
}

fn sorted(members: &HashSet<usize>) -> Vec<usize> {
    let mut positions = members.iter().copied().collect::<Vec<_>>();
    positions.sort_unstable();
    positions
}

/// The graph objective ranks every rung about half as fine as the truth, so where an annotation
/// exists the markers judge the ladder instead.
pub fn pick_rung(ladder: Vec<Partitioning>, judge: &Judge) -> Partitioning {
    let bar = judge.quality.completeness_bar(judge.completeness);
    let passing = |held: &Partitioning| {
        held.cluster_map
            .values()
            .map(|members| judge.score(&sorted(members)))
            .filter(|held| held.completeness >= bar && held.contamination <= TIER_CONTAMINATION)
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
    info!(
        "Markers chose a {} community rung, {value:.0} bins over the bar.",
        chosen.cluster_map.len()
    );
    chosen
}

struct Ranked {
    worth: f64,
    positions: Vec<usize>,
    size_tie: bool,
}

impl PartialEq for Ranked {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Ranked {}

impl PartialOrd for Ranked {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Ranked {
    fn cmp(&self, other: &Self) -> Ordering {
        self.worth
            .total_cmp(&other.worth)
            .then_with(|| match self.size_tie || other.size_tie {
                true => self.positions.len().cmp(&other.positions.len()),
                false => Ordering::Equal,
            })
            .then_with(|| other.positions.cmp(&self.positions))
    }
}

fn candidates(ladder: &[Partitioning]) -> Vec<Vec<usize>> {
    let mut found = ladder
        .iter()
        .flat_map(|held| held.cluster_map.values().map(sorted))
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
    let proposed = candidates(&ladder);
    let mut held = proposed
        .into_iter()
        .filter(|positions| judge.admits(positions))
        .map(|positions| Ranked {
            worth: judge.worth(&positions),
            positions,
            size_tie: judge.size_tie,
        })
        .collect::<BinaryHeap<_>>();

    let mut cluster_map: HashMap<usize, HashSet<usize>> = HashMap::new();
    let mut claimed = HashSet::new();
    while let Some(entry) = held.pop() {
        let left = remaining(&entry.positions, &claimed);
        if left.is_empty() {
            continue;
        }
        if left.len() < entry.positions.len() {
            held.push(Ranked {
                worth: judge.worth(&left),
                positions: left,
                size_tie: judge.size_tie,
            });
            continue;
        }
        claimed.extend(left.iter().copied());
        cluster_map.insert(cluster_map.len(), left.into_iter().collect());
    }

    let outliers = every.difference(&claimed).copied().collect::<HashSet<_>>();
    info!(
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
