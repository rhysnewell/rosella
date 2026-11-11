use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap, HashSet};

use log::debug;
use rayon::prelude::*;

use crate::clustering::clusterer::Partitioning;
use crate::clustering::graph_partition::Partition;
use crate::embedding::features::ContigFeatures;
use crate::markers::ContigMarkers;
use crate::quality::{Quality, Scorer, edge_spread};
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
    pub depth: &'a ContigFeatures<'a>,
    pub markers: &'a ContigMarkers,
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
        self.clears(self.score(positions))
    }

    fn clears(&self, quality: Quality) -> bool {
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
        let taken_contigs = self.mapped(taken);
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
                let found = self.mapped(&offered.found[*rival]);
                let set = self.quality.score(&found).set;
                let worth = self
                    .quality
                    .score_on(&taken_contigs, set)
                    .score(self.bars.worth);
                let spread = edge_spread([
                    (1.0, self.quality.points_on(&found, self.bars.worth, set)),
                    (
                        -1.0,
                        self.quality.points_on(&taken_contigs, self.bars.worth, set),
                    ),
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

// At equal worth the set more binners offer whole goes first, since marker-free contigs leave
// worth level. Under the reporting floor fused fragments outscore each one, so the smaller wins.
struct Entry {
    ranked: Ranked<usize>,
    support: usize,
    below: Option<usize>,
    standing: Standing,
}

// Markers say little about a bin that fails, so among several binners the run's own failing bins
// go ahead of every other binner's, as a run on its own takes them, and then its eukaryotes whole.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Standing {
    Rest,
    HomeEukaryote,
    Home,
    Passing,
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
        let by_worth = || {
            self.ranked
                .worth
                .total_cmp(&other.ranked.worth)
                .then(self.support.cmp(&other.support))
                .then_with(|| self.ranked.cmp(&other.ranked))
        };
        self.standing
            .cmp(&other.standing)
            .then_with(|| match (self.below, other.below) {
                (Some(mine), Some(theirs)) => theirs.cmp(&mine).then_with(by_worth),
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, None) => by_worth(),
            })
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
pub fn combine(ladder: &[Partitioning], voters: &[usize], judge: &Judge) -> Combined {
    let _timer = crate::timing::scope("combine");
    let every = ladder
        .iter()
        .flat_map(|held| held.cluster_map.values().flatten().copied())
        .chain(ladder.iter().flat_map(|held| held.outliers.iter().copied()))
        .collect::<HashSet<_>>();
    let (found, arms): (Vec<_>, Vec<_>) = candidates(ladder).into_iter().unzip();
    let ballots = Ballots::new(&arms, voters);
    let several = ballots.several();
    let mapped = found
        .iter()
        .map(|contigs| judge.mapped(contigs))
        .collect::<Vec<_>>();
    let eukaryotic = judge
        .markers
        .eukaryotic(mapped.iter().map(Vec::as_slice), judge.bars.quality());
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
    // Prokaryote markers rank a eukaryote's pieces by chance, so eukaryotic offers keep their
    // places in the order. The run's own go first by size, since binners shred a eukaryote alike.
    let mut order = offered.worth.clone();
    if several {
        let mut places = (0..offered.found.len())
            .filter(|at| eukaryotic[*at])
            .map(|at| {
                let bases = judge.bases(offered.found[at].iter().copied());
                let home = ballots.home(at);
                let binners = match home {
                    true => 0,
                    false => ballots.agreement(judge, &offered, at),
                };
                (home, binners, bases, Reverse(at))
            })
            .collect::<Vec<_>>();
        places.sort_unstable_by_key(|key| Reverse(*key));
        let mut worths = places
            .iter()
            .map(|(.., Reverse(at))| offered.worth[*at])
            .collect::<Vec<_>>();
        worths.sort_unstable_by(|one, other| other.total_cmp(one));
        for ((.., Reverse(at)), worth) in places.iter().zip(worths) {
            order[*at] = worth;
        }
    }
    let floor = judge.bars.reported().completeness;
    // A eukaryote reads under the reporting floor on prokaryote markers even whole, and another
    // binner's smallest offer under it is a fragment of a bin. The run's own are bins of their own.
    let ranked = |worth: f64, contigs: Vec<usize>, at: usize, support: usize| {
        let quality = judge.score(&contigs);
        let home = ballots.home(at);
        let standing = match (several, eukaryotic[at]) {
            (false, _) => Standing::Passing,
            (true, _) if judge.clears(quality) => Standing::Passing,
            (true, false) if home => Standing::Home,
            (true, true) if home => Standing::HomeEukaryote,
            (true, _) => Standing::Rest,
        };
        let under = match standing {
            Standing::Home => quality.score(judge.bars.worth) < floor,
            _ => !several && !eukaryotic[at] && quality.completeness < floor,
        };
        Entry {
            below: under.then(|| judge.bases(contigs.iter().copied())),
            ranked: Ranked {
                worth,
                contigs,
                extra: at,
            },
            support,
            standing,
        }
    };
    // A eukaryote's place was not set by its markers, so what is left of one keeps the place.
    let entry = |contigs: Vec<usize>, at: usize, support: usize, place: f64| {
        let worth = match several && eukaryotic[at] {
            true => place,
            false => judge.worth(&contigs),
        };
        ranked(worth, contigs, at, support)
    };
    let mut held = offered
        .found
        .iter()
        .zip(&order)
        .enumerate()
        .map(|(at, (contigs, worth))| ranked(*worth, contigs.clone(), at, ballots.support(at)))
        .collect::<BinaryHeap<_>>();

    let mut cluster_map = BTreeMap::new();
    let mut sources = Vec::new();
    let mut claimed = HashSet::new();
    let mut free = offered
        .found
        .iter()
        .map(|contigs| judge.bases(contigs.iter().copied()))
        .collect::<Vec<_>>();
    while let Some(Entry {
        ranked: popped,
        support,
        below,
        ..
    }) = held.pop()
    {
        let left = remaining(&popped.contigs, &claimed);
        if left.is_empty() {
            continue;
        }
        if left.len() < popped.contigs.len() {
            held.push(entry(left, popped.extra, support, popped.worth));
            continue;
        }
        let widened = match below {
            Some(_) => None,
            None => judge.wider(
                &left,
                popped.worth,
                &offered,
                &claimed,
                |rival| ballots.support(rival),
                support,
            ),
        };
        let (at, taken) = match widened {
            Some(rival) => (rival, offered.found[rival].clone()),
            None => (popped.extra, left),
        };
        let passes = several && judge.passes(&taken);
        let pieces = match several && !passes && !eukaryotic[at] {
            true => ballots.pieces(judge, &offered, &taken, &claimed, &free, |rival| {
                rival != at && !eukaryotic[rival]
            }),
            false => Vec::new(),
        };
        let taken = match (passes, pieces.is_empty()) {
            (true, _) => vec![(
                at,
                ballots
                    .resolve(judge, &offered, at, &taken, &claimed)
                    .unwrap_or(taken),
            )],
            (false, true) => vec![(at, taken)],
            (false, false) => pieces,
        };
        for (at, taken) in taken {
            for contig in taken.iter().filter(|contig| claimed.insert(**contig)) {
                let bases = judge.bases([*contig]);
                for holder in &offered.holders[contig] {
                    free[*holder] -= bases;
                }
            }
            sources.push((arms[at].clone(), taken.len() == offered.found[at].len()));
            cluster_map.insert(cluster_map.len(), taken);
        }
        let rest = remaining(&popped.contigs, &claimed);
        if !rest.is_empty() {
            held.push(entry(rest, popped.extra, support, popped.worth));
        }
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
