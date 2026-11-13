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
use crate::refine::ranking::{Ranked, kept, sorted};
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
        backed: impl Fn(usize) -> bool,
    ) -> Option<usize> {
        if worth <= 0.0 {
            return None;
        }
        let size = self.bases(taken.iter().copied());
        let held = taken.iter().copied().collect::<HashSet<_>>();
        taken
            .iter()
            .flat_map(|at| offered.holders[at].iter().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .filter(|rival| offered.worth[*rival] > 0.0 && backed(*rival))
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
            .filter(|(rival, _)| self.level(taken, offered, *rival))
            .max_by_key(|(rival, bases)| (*bases, Reverse(*rival)))
            .map(|(rival, _)| rival)
    }

    fn level(&self, taken: &[usize], offered: &Offered, rival: usize) -> bool {
        let taken = self.mapped(taken);
        let found = self.mapped(&offered.found[rival]);
        let set = self.quality.score(&found).set;
        let worth = self.quality.score_on(&taken, set).score(self.bars.worth);
        let spread = edge_spread([
            (1.0, self.quality.points_on(&found, self.bars.worth, set)),
            (-1.0, self.quality.points_on(&taken, self.bars.worth, set)),
        ]);
        offered.worth[rival].powi(2) - worth.powi(2) >= -spread
    }
}

pub(crate) struct Offered {
    pub(crate) found: Vec<Vec<usize>>,
    worth: Vec<f64>,
    pub(crate) holders: HashMap<usize, Vec<usize>>,
}

pub struct Combined {
    pub partitioning: Partitioning,
    pub offers: Vec<Offer>,
    pub picks: Vec<Pick>,
}

pub struct Offer {
    pub arms: Vec<usize>,
    pub first: Option<usize>,
    pub contigs: usize,
    pub quality: Quality,
    pub worth: f64,
    pub support: usize,
    pub standing: Standing,
}

pub struct Pick {
    pub popped: usize,
    pub taken: usize,
    pub worth: f64,
    pub contigs: usize,
}

// Marker-free contigs leave worth level, so a tie goes to the set more binners offer, then the
// run's own. Under the floor fused fragments outscore each one, so the smaller wins.
struct Entry {
    ranked: Ranked<usize>,
    support: usize,
    below: Option<usize>,
    standing: Standing,
    home: bool,
    version: usize,
}

fn unclaimed(contigs: &[usize], taken: &[bool]) -> Vec<usize> {
    kept(contigs, |at| !taken[*at])
}

fn requeue(held: &mut BinaryHeap<Entry>, version: &mut [usize], mut fresh: Entry) {
    let at = fresh.ranked.extra;
    version[at] += 1;
    fresh.version = version[at];
    held.push(fresh);
}

// Markers say little about a bin that fails, so among several binners the run's own failing bins
// go ahead of every other binner's, as a run on its own takes them, and then its eukaryotes whole.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Standing {
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
                .then(self.home.cmp(&other.home))
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
    let scores = found
        .par_iter()
        .map(|contigs| judge.score(contigs))
        .collect::<Vec<_>>();
    let worth = scores
        .iter()
        .map(|quality| quality.score(judge.bars.worth))
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
    let ranked = |worth: f64, quality: Quality, contigs: Vec<usize>, at: usize, support: usize| {
        let home = ballots.home(at);
        let standing = match (several, eukaryotic[at]) {
            (false, _) => Standing::Passing,
            (true, _) if judge.clears(quality) => Standing::Passing,
            (true, false) if home => Standing::Home,
            (true, true) if home => Standing::HomeEukaryote,
            (true, _) => Standing::Rest,
        };
        // A bag of small genomes reads over the floor on completeness while its contamination
        // sinks its worth, so the floor is read on worth or the bag pops before their own bins.
        let judged = standing == Standing::Home || (!several && !eukaryotic[at]);
        let under = judged && quality.score(judge.bars.worth) < floor;
        Entry {
            below: under.then(|| judge.bases(contigs.iter().copied())),
            ranked: Ranked {
                worth,
                contigs,
                extra: at,
            },
            support,
            standing,
            home: several && home,
            version: 0,
        }
    };
    // A eukaryote's place was not set by its markers, so what is left of one keeps the place.
    let read = |contigs: Vec<usize>, at: usize, support: usize, place: f64, quality: Quality| {
        let worth = match several && eukaryotic[at] {
            true => place,
            false => quality.score(judge.bars.worth),
        };
        ranked(worth, quality, contigs, at, support)
    };
    let entry = |contigs: Vec<usize>, at: usize, support: usize, place: f64| {
        let quality = judge.score(&contigs);
        (quality, read(contigs, at, support, place, quality))
    };
    let entries = offered
        .found
        .iter()
        .zip(&order)
        .zip(&scores)
        .enumerate()
        .map(|(at, ((contigs, worth), quality))| {
            ranked(*worth, *quality, contigs.clone(), at, ballots.support(at))
        })
        .collect::<Vec<_>>();
    let offers = entries
        .iter()
        .zip(&scores)
        .map(|(entry, quality)| Offer {
            arms: arms[entry.ranked.extra].clone(),
            first: entry.ranked.contigs.first().copied(),
            contigs: entry.ranked.contigs.len(),
            quality: *quality,
            worth: entry.ranked.worth,
            support: entry.support,
            standing: entry.standing,
        })
        .collect();
    let mut held = entries.into_iter().collect::<BinaryHeap<_>>();
    let mut version = vec![0; offered.found.len()];
    let mut now = scores.clone();
    let mut taken_at = vec![false; judge.contigs.len()];
    let stale = std::env::var_os("ROSELLA_EUKARYOTE_STALE").is_some();

    let mut cluster_map = BTreeMap::new();
    let mut picks = Vec::new();
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
        version: seen,
        ..
    }) = held.pop()
    {
        if seen != version[popped.extra] {
            continue;
        }
        let left = unclaimed(&popped.contigs, &taken_at);
        if left.is_empty() {
            continue;
        }
        if left.len() < popped.contigs.len() {
            let (quality, fresh) = entry(left, popped.extra, support, popped.worth);
            now[popped.extra] = quality;
            requeue(&mut held, &mut version, fresh);
            continue;
        }
        let widened = match below {
            // Only the run's own bin was cut on this run's depth, so it may widen a winner that
            // more binners offer.
            None => judge.wider(&left, popped.worth, &offered, &claimed, |rival| {
                (several && ballots.home(rival)) || ballots.support(rival) >= support
            }),
            Some(_) => None,
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
        let mut touched = BTreeMap::<usize, Vec<usize>>::new();
        for (at, taken) in taken {
            for contig in taken.iter().filter(|contig| claimed.insert(**contig)) {
                taken_at[*contig] = true;
                let bases = judge.bases([*contig]);
                for holder in &offered.holders[contig] {
                    free[*holder] -= bases;
                    touched.entry(*holder).or_default().push(*contig);
                }
            }
            picks.push(Pick {
                popped: popped.extra,
                taken: at,
                worth: popped.worth,
                contigs: taken.len(),
            });
            cluster_map.insert(cluster_map.len(), taken);
        }
        let rest = unclaimed(&popped.contigs, &taken_at);
        if !rest.is_empty() {
            let (quality, fresh) = entry(rest, popped.extra, support, popped.worth);
            now[popped.extra] = quality;
            requeue(&mut held, &mut version, fresh);
        }
        // A rival that lost another genome's contigs can now read clean over the floor, so it is
        // scored again here. A eukaryote mix that sheds prokaryotes reads clean too. It waits.
        let rescored = touched
            .into_iter()
            .filter(|(rival, _)| *rival != popped.extra && !(stale && eukaryotic[*rival]))
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|(rival, lost)| {
                let left = unclaimed(&offered.found[rival], &taken_at);
                let fresh = (!left.is_empty()).then(|| {
                    let quality = match judge.markers.reads_without(
                        &judge.mapped(&lost),
                        free[rival],
                        now[rival],
                    ) {
                        true => now[rival],
                        false => judge.score(&left),
                    };
                    let support = ballots.support(rival);
                    (quality, read(left, rival, support, order[rival], quality))
                });
                (rival, fresh)
            })
            .collect::<Vec<_>>();
        for (rival, fresh) in rescored {
            match fresh {
                None => version[rival] += 1,
                Some((quality, fresh)) => {
                    now[rival] = quality;
                    requeue(&mut held, &mut version, fresh);
                }
            }
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
        offers,
        picks,
    }
}
