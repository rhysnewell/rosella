use std::cmp::Reverse;
use std::collections::{BTreeMap, HashMap, HashSet};

use log::debug;

use crate::quality::Quality;
use crate::recover::agree::{Agreement, Sums};
use crate::recover::fold::{HOST, Role};
use crate::recover::join_report::Report;
use crate::recover::ladder::{Judge, Offered};
use crate::refine::ranking::sorted;

pub type Group = Vec<Vec<usize>>;

// Under the floor the smaller candidate is taken first, so a genome cut on depth noise leaves as
// pieces. A piece that agrees with the rest of a larger candidate holding it rejoins that one.
pub(crate) struct Joins<'a> {
    agreement: Agreement<'a>,
    sums: HashMap<usize, Sums>,
    gains: bool,
    clean_only: bool,
    joined: Vec<(usize, usize)>,
    taken: Vec<(usize, usize)>,
    report: Option<Report>,
}

impl<'a> Joins<'a> {
    pub(crate) fn of(
        judge: &Judge<'a>,
        offered: &Offered,
        scores: &[Quality],
        eukaryotic: &[bool],
    ) -> Option<Self> {
        let mut hosts = offered
            .found
            .iter()
            .zip(scores.iter().zip(eukaryotic))
            .filter(|(_, (quality, eukaryotic))| Role::of(**quality, **eukaryotic) == Role::Host)
            .map(|(contigs, _)| judge.mapped(contigs))
            .collect::<Vec<_>>();
        hosts.sort_by_key(|members| {
            let bases = members
                .iter()
                .map(|contig| judge.lengths[*contig])
                .sum::<usize>();
            (Reverse(bases), members.clone())
        });
        let mut seen = HashSet::new();
        hosts.retain(|members| {
            let fresh = members.iter().all(|contig| !seen.contains(contig));
            if fresh {
                seen.extend(members.iter().copied());
            }
            fresh
        });
        let agreement = Agreement::of(
            &hosts,
            judge.depth.depths(),
            judge.depth.composition(),
            judge.lengths,
        )?;
        Some(Self {
            agreement,
            sums: HashMap::new(),
            gains: std::env::var("ROSELLA_UNDER_JOIN").as_deref() == Ok("gain"),
            clean_only: std::env::var("ROSELLA_UNDER_JOIN").as_deref() == Ok("clean"),
            joined: Vec::new(),
            taken: Vec::new(),
            report: Report::open(),
        })
    }

    pub(crate) fn target(
        &mut self,
        judge: &Judge,
        offered: &Offered,
        piece: &[usize],
        at: usize,
        eligible: impl Fn(usize) -> bool,
    ) -> Option<usize> {
        let size = judge.bases(piece.iter().copied());
        let mut shared = HashMap::<usize, Vec<usize>>::new();
        for contig in piece {
            for holder in offered.holders[contig]
                .iter()
                .filter(|holder| **holder != at)
            {
                shared.entry(*holder).or_default().push(*contig);
            }
        }
        let mut rivals = shared
            .into_iter()
            .filter(|(rival, _)| eligible(*rival))
            .map(|(rival, shared)| {
                let bases = judge.bases(offered.found[rival].iter().copied());
                (bases, rival, shared)
            })
            .filter(|(bases, _, shared)| {
                *bases > size && 2 * judge.bases(shared.iter().copied()) >= size
            })
            .collect::<Vec<_>>();
        rivals.sort_unstable_by_key(|(bases, rival, _)| (*bases, *rival));
        let mapped = judge.mapped(piece);
        let mine = self.agreement.sums(mapped.iter().copied());
        if let Some(report) = &mut self.report {
            report.piece(&mapped);
        }
        let held = mapped.iter().collect::<HashSet<_>>();
        let mut chosen = None;
        for (bases, rival, shared) in rivals {
            let found = judge.mapped(&offered.found[rival]);
            let whole = self
                .sums
                .entry(rival)
                .or_insert_with(|| self.agreement.sums(found.iter().copied()));
            let rest = whole.without(&self.agreement.sums(judge.mapped(&shared)));
            let agree = self.agreement.agree(&mine, &rest);
            let union = (agree || self.report.is_some()).then(|| {
                let held = shared.iter().collect::<HashSet<_>>();
                let outside = piece.iter().filter(|at| !held.contains(at)).copied();
                found
                    .iter()
                    .copied()
                    .chain(judge.mapped(&outside.collect::<Vec<_>>()))
                    .collect::<Vec<_>>()
            });
            let joined = union.as_deref().map(|union| judge.quality.score(union));
            let rest_contigs = union
                .as_deref()
                .unwrap_or_default()
                .iter()
                .copied()
                .filter(|contig| !held.contains(contig))
                .collect::<Vec<_>>();
            let host = joined.is_some_and(|quality| match self.clean_only {
                true => quality.contamination <= judge.bars.tier(),
                false => quality.clears(HOST),
            });
            let rest_quality = (host && self.gains || self.report.is_some())
                .then(|| judge.quality.score(&rest_contigs));
            let gains = !self.gains
                || rest_quality
                    .zip(joined)
                    .is_some_and(|(rest, joined)| joined.completeness > rest.completeness);
            let pick = chosen.is_none() && agree && host && gains;
            if pick {
                chosen = Some(rival);
            }
            if let Some(report) = &mut self.report {
                let statistics = self.agreement.statistics(&mine, &rest);
                let p = judge.quality.score(&mapped);
                let r = rest_quality.unwrap_or(p);
                let u = joined.unwrap_or(p);
                let fields = [
                    size as f64,
                    bases as f64,
                    judge.bases(shared.iter().copied()) as f64,
                    statistics.depth.0,
                    statistics.depth.1 as f64,
                    statistics.composition.0,
                    statistics.composition.1 as f64,
                    p.completeness,
                    p.contamination,
                    r.completeness,
                    r.contamination,
                    u.completeness,
                    u.contamination,
                    agree as u8 as f64,
                    pick as u8 as f64,
                ]
                .map(|value| format!("{value:.4}"));
                report.test(rival, &found, &fields);
            } else if chosen.is_some() {
                break;
            }
        }
        chosen
    }

    pub(crate) fn record(&mut self, cluster: usize, at: usize, target: Option<usize>) {
        self.taken.push((cluster, at));
        if let Some(target) = target {
            self.joined.push((cluster, target));
        }
    }

    pub(crate) fn groups(&self, cluster_map: &BTreeMap<usize, Vec<usize>>) -> Vec<Group> {
        let mut parent = (0..cluster_map.len()).collect::<Vec<_>>();
        let mut first = HashMap::<usize, usize>::new();
        for (cluster, target) in &self.joined {
            let head = *first.entry(*target).or_insert(*cluster);
            union(&mut parent, head, *cluster);
        }
        for (cluster, at) in &self.taken {
            if let Some(head) = first.get(at) {
                union(&mut parent, *head, *cluster);
            }
        }
        let mut groups = BTreeMap::<usize, Group>::new();
        for (cluster, members) in cluster_map {
            let head = root(&mut parent, *cluster);
            groups
                .entry(head)
                .or_default()
                .push(sorted(members.iter().copied()));
        }
        let groups = groups
            .into_values()
            .filter(|group| group.len() > 1)
            .collect::<Vec<_>>();
        debug!(
            "Under-floor joins: {} pieces onto {} candidates in {} groups",
            self.joined.len(),
            first.len(),
            groups.len()
        );
        groups
    }
}

// The weight walk and the pass contest judge each pass on its own bins, so a pass's joins land
// only where every bin they join came through to the end.
pub fn joined(
    cluster_map: BTreeMap<usize, Vec<usize>>,
    groups: &[Group],
) -> BTreeMap<usize, Vec<usize>> {
    let bins = cluster_map.into_values().collect::<Vec<_>>();
    let at = bins
        .iter()
        .enumerate()
        .map(|(at, members)| (sorted(members.iter().copied()), at))
        .collect::<HashMap<_, _>>();
    let mut parent = (0..bins.len()).collect::<Vec<_>>();
    for group in groups {
        let found = group
            .iter()
            .map(|members| at.get(members).copied())
            .collect::<Option<Vec<_>>>();
        if let Some(found) = found {
            for pair in found.windows(2) {
                union(&mut parent, pair[0], pair[1]);
            }
        }
    }
    let mut merged = BTreeMap::<usize, Vec<usize>>::new();
    for (at, members) in bins.into_iter().enumerate() {
        let head = root(&mut parent, at);
        merged.entry(head).or_default().extend(members);
    }
    merged
        .into_values()
        .enumerate()
        .map(|(id, members)| (id, sorted(members)))
        .collect()
}

fn root(parent: &mut [usize], mut at: usize) -> usize {
    while parent[at] != at {
        parent[at] = parent[parent[at]];
        at = parent[at];
    }
    at
}

fn union(parent: &mut [usize], one: usize, other: usize) {
    let (one, other) = (root(parent, one), root(parent, other));
    let (low, high) = (one.min(other), one.max(other));
    parent[high] = low;
}
