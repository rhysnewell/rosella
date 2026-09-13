use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};

use log::{debug, info};

use crate::clustering::clusterer::Partitioning;
use crate::clustering::graph_partition::Partition;
use crate::quality::{Quality, Scorer, Worth};
use crate::refine::select::remaining;

pub const COMBINE_SOURCE_NAMES: [&str; 2] = ["ladder", "arms"];

/// The fine rungs are where the fragments that cut a whole genome come from, so combining can be
/// offered the whole ladder or only the objective's pick from each partition arm.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum CombineSource {
    #[default]
    Ladder,
    Arms,
}

impl CombineSource {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "ladder" => Some(Self::Ladder),
            "arms" => Some(Self::Arms),
            _ => None,
        }
    }
}

/// The ladder arrives sorted by the graph objective, so the first entry an arm contributes is
/// that arm's best.
pub fn best_per_arm(ladder: Vec<Partitioning>) -> Vec<Partitioning> {
    let mut seen = Vec::new();
    ladder
        .into_iter()
        .filter(|held| match seen.contains(&held.arm) {
            true => false,
            false => {
                seen.push(held.arm);
                true
            }
        })
        .collect()
}

pub const RUNG_STATISTIC_NAMES: [&str; 4] = ["pass50", "pass80", "pass90", "worth"];

/// What ranks one rung of the ladder against another. Counting communities over a completeness
/// bar turns over near the true genome count; summing worth climbs with community count.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum RungStatistic {
    #[default]
    Pass50,
    Pass80,
    Pass90,
    Worth,
}

impl RungStatistic {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "pass50" => Some(Self::Pass50),
            "pass80" => Some(Self::Pass80),
            "pass90" => Some(Self::Pass90),
            "worth" => Some(Self::Worth),
            _ => None,
        }
    }
}

pub struct Judge<'a> {
    pub quality: &'a dyn Scorer,
    pub contigs: &'a [usize],
    pub worth: Worth,
    pub rung_statistic: RungStatistic,
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
}

fn sorted(members: &HashSet<usize>) -> Vec<usize> {
    let mut positions = members.iter().copied().collect::<Vec<_>>();
    positions.sort_unstable();
    positions
}

/// The graph objective ranks every rung about half as fine as the truth, so where an annotation
/// exists the markers judge the ladder instead.
pub fn pick_rung(ladder: Vec<Partitioning>, judge: &Judge) -> Partitioning {
    let worth = |held: &Partitioning| {
        let scored = held
            .cluster_map
            .values()
            .map(|members| judge.score(&sorted(members)))
            .collect::<Vec<_>>();
        let sum = scored
            .iter()
            .map(|held| held.score(judge.worth).max(0.0))
            .sum::<f64>();
        let clean = |bar: f64| {
            scored
                .iter()
                .filter(|held| held.completeness >= bar && held.contamination <= 10.0)
                .count()
        };
        (sum, clean(50.0), clean(80.0), clean(90.0))
    };
    let mut scored = ladder
        .into_iter()
        .map(|held| {
            let (sum, medium, upper, high) = worth(&held);
            debug!(
                "rung {} communities sum {sum:.1} pass50 {medium} pass80 {upper} pass90 {high}",
                held.cluster_map.len()
            );
            let ranked = match judge.rung_statistic {
                RungStatistic::Pass50 => medium as f64,
                RungStatistic::Pass80 => upper as f64,
                RungStatistic::Pass90 => high as f64,
                RungStatistic::Worth => sum,
            };
            (ranked, held)
        })
        .collect::<Vec<_>>();
    scored.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(Ordering::Equal));
    let (value, chosen) = scored.swap_remove(0);
    info!(
        "Markers chose a {} community rung, {:?} {value:.1}.",
        chosen.cluster_map.len(),
        judge.rung_statistic
    );
    chosen
}

struct Ranked {
    worth: f64,
    positions: Vec<usize>,
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
        .map(|positions| Ranked {
            worth: judge.worth(&positions),
            positions,
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
        score: f64::NAN,
        arm: Partition::Both,
    }
}
