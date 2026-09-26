use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};

use log::debug;

use crate::clustering::clusterer::Partitioning;
use crate::clustering::graph_partition::Partition;
use crate::quality::Quality;
use crate::recover::combine_report::CombineReport;
use crate::recover::ladder::{Judge, candidates};
use crate::refine::select::remaining;

const TIERS: [f64; 6] = [10.0, 20.0, 30.0, 40.0, 50.0, 100.0];

/// F1 is multiplicative, so completeness cannot buy contamination the way it can under worth.
fn merit(held: &Quality) -> f64 {
    let recall = held.completeness / 100.0;
    let precision = (1.0 - held.contamination / 100.0).max(0.0);
    match recall + precision > 0.0 {
        true => 2.0 * recall * precision / (recall + precision),
        false => 0.0,
    }
}

struct Peeled {
    merit: f64,
    bp: usize,
    contigs: Vec<usize>,
}

impl PartialEq for Peeled {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Peeled {}

impl PartialOrd for Peeled {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Peeled {
    fn cmp(&self, other: &Self) -> Ordering {
        self.merit
            .total_cmp(&other.merit)
            .then_with(|| other.bp.cmp(&self.bp))
            .then_with(|| other.contigs.cmp(&self.contigs))
    }
}

fn bases(positions: &[usize], judge: &Judge<'_>, lengths: &[usize]) -> usize {
    positions
        .iter()
        .filter_map(|at| judge.contigs.get(*at))
        .filter_map(|contig| lengths.get(*contig))
        .sum()
}

fn admit(
    contigs: Vec<usize>,
    tier: f64,
    judge: &Judge<'_>,
    lengths: &[usize],
    heap: &mut BinaryHeap<Peeled>,
    deferred: &mut Vec<Vec<usize>>,
) {
    let held = judge.score(&contigs);
    match held.contamination <= tier {
        true => heap.push(Peeled {
            merit: merit(&held),
            bp: bases(&contigs, judge, lengths),
            contigs,
        }),
        false => deferred.push(contigs),
    }
}

/// Worth ranks a 95 complete, 10 contaminated candidate over a 60 complete, clean one, and then
/// the first one claims the second one's contigs. Draining the ladder in contamination tiers is
/// a claim order rather than a bar: a pure candidate takes what it owns before a dirtier one is
/// offered the same contigs at all.
pub fn peel(
    ladder: &[Partitioning],
    judge: &Judge<'_>,
    lengths: &[usize],
    report: Option<&CombineReport<'_>>,
) -> Partitioning {
    let _timer = crate::timing::scope("peel");
    let every = ladder
        .iter()
        .flat_map(|held| held.cluster_map.values().flatten().copied())
        .chain(ladder.iter().flat_map(|held| held.outliers.iter().copied()))
        .collect::<HashSet<_>>();

    let mut waiting = candidates(ladder);
    let mut cluster_map: HashMap<usize, HashSet<usize>> = HashMap::new();
    let mut claimed: HashSet<usize> = HashSet::new();
    let mut seen = 0;

    for tier in TIERS.into_iter().chain(std::iter::once(f64::INFINITY)) {
        let mut heap = BinaryHeap::new();
        let mut deferred = Vec::new();
        for contigs in waiting.drain(..) {
            let left = remaining(&contigs, &claimed);
            if left.is_empty() {
                continue;
            }
            admit(left, tier, judge, lengths, &mut heap, &mut deferred);
        }
        while let Some(entry) = heap.pop() {
            let left = remaining(&entry.contigs, &claimed);
            if let Some(report) = report {
                let verdict = match (left.is_empty(), left.len() < entry.contigs.len()) {
                    (true, _) => "spent",
                    (_, true) => "trimmed",
                    _ => "taken",
                };
                report.row(seen, entry.merit, verdict, &judge.mapped(&entry.contigs));
                seen += 1;
            }
            if left.is_empty() {
                continue;
            }
            if left.len() < entry.contigs.len() {
                admit(left, tier, judge, lengths, &mut heap, &mut deferred);
                continue;
            }
            claimed.extend(left.iter().copied());
            cluster_map.insert(cluster_map.len(), left.into_iter().collect());
        }
        waiting = deferred;
    }

    let outliers = every.difference(&claimed).copied().collect::<HashSet<_>>();
    debug!(
        "Markers peeled {} rungs into {} bins, {} contigs unclaimed.",
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
