use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap, HashSet};

use log::debug;
use rayon::prelude::*;

use crate::clustering::clusterer::Partitioning;
use crate::clustering::graph_partition::Partition;
use crate::embedding::features::ContigFeatures;
use crate::quality::{Quality, Scorer, edge_spread};
use crate::recover::combine_report::CombineReport;
use crate::recover::consensus::Ballots;
use crate::refine::ranking::{Ranked, remaining, sorted};
use crate::refine::rung::Bars;

// Two neighbouring rungs often clear the bar within a few bins of each other, and each keeps whole
// genomes the other cuts, so both go to the arbiter. A third only adds fused candidates.
const KEPT_RUNGS: usize = 2;

pub fn top_rungs_per_arm(ladder: Vec<Partitioning>, judge: &Judge) -> Vec<Partitioning> {
    let mut arms: Vec<((Partition, u64), Vec<Partitioning>)> = Vec::new();
    for held in ladder {
        let key = (held.arm, held.seed);
        match arms.iter_mut().find(|(seen, _)| *seen == key) {
            Some((_, entries)) => entries.push(held),
            None => arms.push((key, vec![held])),
        }
    }
    arms.into_iter()
        .flat_map(|(_, entries)| top_rungs(entries, judge))
        .collect()
}

pub struct Judge<'a> {
    pub quality: &'a dyn Scorer,
    pub contigs: &'a [usize],
    pub lengths: &'a [usize],
    pub bars: Bars,
    pub depth: Option<&'a ContigFeatures<'a>>,
}

impl Judge<'_> {
    pub(crate) fn score(&self, positions: &[usize]) -> Quality {
        self.quality.score(&self.mapped(positions))
    }

    fn worth(&self, positions: &[usize]) -> f64 {
        self.score(positions).score(self.bars.worth)
    }

    pub(crate) fn mapped(&self, positions: &[usize]) -> Vec<usize> {
        positions.iter().map(|at| self.contigs[*at]).collect()
    }

    pub(crate) fn passes(&self, positions: &[usize]) -> bool {
        let quality = self.score(positions);
        let floor = self.bars.reported().completeness;
        quality.completeness >= floor
            && quality.contamination <= self.bars.tier()
            && quality.score(self.bars.worth) >= floor
    }

    pub(crate) fn bases(&self, positions: impl IntoIterator<Item = usize>) -> usize {
        positions
            .into_iter()
            .map(|at| self.lengths[self.contigs[at]])
            .sum()
    }

    fn points(&self, positions: &[usize]) -> Vec<(usize, f64)> {
        self.quality
            .points(&self.mapped(positions), self.bars.worth)
    }

    // The markers cannot see a genome's own marker-free contigs, so a larger candidate that holds
    // most of the winner and reads level with it inside one standard error keeps them.
    fn wider(
        &self,
        taken: &[usize],
        worth: f64,
        offered: &Offered,
        claimed: &HashSet<usize>,
        support: impl Fn(usize) -> usize,
        least: usize,
    ) -> Option<usize> {
        if worth <= 0.0 {
            return None;
        }
        let size = self.bases(taken.iter().copied());
        let held = taken.iter().copied().collect::<HashSet<_>>();
        let points = self.points(taken);
        taken
            .iter()
            .flat_map(|at| offered.holders[at].iter().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .filter(|rival| offered.worth[*rival] > 0.0 && support(*rival) >= least)
            .filter(|rival| offered.found[*rival].iter().all(|at| !claimed.contains(at)))
            .map(|rival| (rival, self.bases(offered.found[rival].iter().copied())))
            .filter(|(rival, bases)| {
                *bases > size
                    && 2 * self.bases(
                        offered.found[*rival]
                            .iter()
                            .copied()
                            .filter(|at| held.contains(at)),
                    ) >= size
            })
            .filter(|(rival, _)| {
                let spread = edge_spread([
                    (1.0, self.points(&offered.found[*rival])),
                    (-1.0, points.clone()),
                ]);
                offered.worth[*rival].powi(2) - worth.powi(2) >= -spread
            })
            .max_by_key(|(rival, bases)| (*bases, Reverse(*rival)))
            .map(|(rival, _)| rival)
    }
}

pub(crate) struct Offered {
    pub(crate) found: Vec<Vec<usize>>,
    worth: Vec<f64>,
    pub(crate) holders: HashMap<usize, Vec<usize>>,
}

pub struct Combined {
    pub partitioning: Partitioning,
    pub sources: Vec<(Vec<usize>, bool)>,
}

// Marker-free contigs leave worth level, so at equal worth the set more binners offer whole
// goes first, rather than one binner's copy carrying contigs no other binner keeps there.
struct Entry {
    ranked: Ranked<usize>,
    support: usize,
}

impl PartialEq for Entry {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}

impl Eq for Entry {}

impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Entry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.ranked
            .worth
            .total_cmp(&other.ranked.worth)
            .then(self.support.cmp(&other.support))
            .then_with(|| self.ranked.cmp(&other.ranked))
    }
}

impl Entry {
    fn new(worth: f64, contigs: Vec<usize>, at: usize, support: usize) -> Self {
        Self {
            ranked: Ranked {
                worth,
                contigs,
                extra: at,
            },
            support,
        }
    }
}

// The ladder arrives in codelength order and the sort is stable, so codelength breaks a tie.
pub fn top_rungs(ladder: Vec<Partitioning>, judge: &Judge) -> Vec<Partitioning> {
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
                "Rung: {} communities ({count} over the bar)",
                held.cluster_map.len()
            );
            (count, held)
        })
        .collect::<Vec<_>>();
    scored.sort_by_key(|(count, _)| Reverse(*count));
    scored
        .into_iter()
        .take(KEPT_RUNGS)
        .map(|(_, held)| held)
        .collect()
}

pub(crate) fn candidates(ladder: &[Partitioning]) -> Vec<(Vec<usize>, Vec<usize>)> {
    let mut found = BTreeMap::<Vec<usize>, Vec<usize>>::new();
    for (arm, held) in ladder.iter().enumerate() {
        for members in held.cluster_map.values() {
            let arms = found.entry(members.clone()).or_default();
            if arms.last() != Some(&arm) {
                arms.push(arm);
            }
        }
    }
    found.into_iter().collect()
}

// Choosing one rung whole is worth almost nothing against choosing the best of both arms, so the
// bins are arbitrated one at a time instead and a rung contributes only the ones that win.
pub fn combine_traced(
    ladder: &[Partitioning],
    voters: &[usize],
    judge: &Judge,
    report: Option<&CombineReport<'_>>,
) -> Combined {
    let _timer = crate::timing::scope("combine");
    let every = ladder
        .iter()
        .flat_map(|held| held.cluster_map.values().flatten().copied())
        .chain(ladder.iter().flat_map(|held| held.outliers.iter().copied()))
        .collect::<HashSet<_>>();
    let (found, arms): (Vec<_>, Vec<_>) = candidates(ladder).into_iter().unzip();
    let ballots = Ballots::new(&arms, voters);
    let mut holders = HashMap::<usize, Vec<usize>>::new();
    for (at, contigs) in found.iter().enumerate() {
        for contig in contigs {
            holders.entry(*contig).or_default().push(at);
        }
    }
    let worth = found
        .par_iter()
        .map(|contigs| judge.worth(contigs))
        .collect::<Vec<_>>();
    let offered = Offered {
        found,
        worth,
        holders,
    };
    let mut held = offered
        .found
        .iter()
        .zip(&offered.worth)
        .enumerate()
        .map(|(at, (contigs, worth))| Entry::new(*worth, contigs.clone(), at, ballots.support(at)))
        .collect::<BinaryHeap<_>>();

    let mut cluster_map = BTreeMap::new();
    let mut sources = Vec::new();
    let mut claimed = HashSet::new();
    let mut seen = 0;
    let mut row = |worth: f64, verdict: &str, contigs: &[usize]| {
        if let Some(report) = report {
            report.row(seen, worth, verdict, &judge.mapped(contigs));
            seen += 1;
        }
    };
    while let Some(Entry {
        ranked: entry,
        support,
    }) = held.pop()
    {
        let left = remaining(&entry.contigs, &claimed);
        if left.is_empty() {
            row(entry.worth, "spent", &entry.contigs);
            continue;
        }
        if left.len() < entry.contigs.len() {
            row(entry.worth, "trimmed", &entry.contigs);
            held.push(Entry::new(judge.worth(&left), left, entry.extra, support));
            continue;
        }
        let widened = judge.wider(
            &left,
            entry.worth,
            &offered,
            &claimed,
            |rival| ballots.support(rival),
            support,
        );
        let (at, taken) = match widened {
            Some(rival) => {
                row(entry.worth, "widened", &entry.contigs);
                row(offered.worth[rival], "taken", &offered.found[rival]);
                (rival, offered.found[rival].clone())
            }
            None => {
                row(entry.worth, "taken", &entry.contigs);
                (entry.extra, left)
            }
        };
        let whole = taken.len() == offered.found[at].len();
        let taken = match ballots.resolve(judge, &offered, at, &taken, &claimed) {
            Some(voted) => {
                row(judge.worth(&voted), "voted", &voted);
                voted
            }
            None => taken,
        };
        claimed.extend(taken.iter().copied());
        let rest = remaining(&entry.contigs, &claimed);
        if !rest.is_empty() {
            held.push(Entry::new(judge.worth(&rest), rest, entry.extra, support));
        }
        sources.push((arms[at].clone(), whole));
        cluster_map.insert(cluster_map.len(), taken);
    }

    let outliers = sorted(every.difference(&claimed).copied());
    debug!(
        "Combined: {} rungs into {} bins ({} contigs unclaimed)",
        ladder.len(),
        cluster_map.len(),
        outliers.len()
    );
    Combined {
        partitioning: Partitioning {
            cluster_map,
            outliers,
            score: None,
            arm: Partition::Both,
            seed: 0,
        },
        sources,
    }
}
