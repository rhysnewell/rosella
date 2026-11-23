use std::cmp::Reverse;
use std::collections::{BTreeMap, HashMap};

use log::debug;

use crate::quality::Quality;
use crate::recover::agree::{Agreement, Sums};
use crate::recover::fold::Role;
use crate::recover::ladder::{Judge, Offered};

// Under the floor the smaller candidate is taken first, so a genome cut on depth noise leaves as
// pieces. A piece that agrees with the rest of a larger candidate holding it rejoins that one.
pub(crate) struct Joins<'a> {
    agreement: Agreement<'a>,
    sums: HashMap<usize, Sums>,
    joined: Vec<(usize, usize)>,
    taken: Vec<(usize, usize)>,
}

impl<'a> Joins<'a> {
    pub(crate) fn of(
        judge: &Judge<'a>,
        offered: &Offered,
        scores: &[Quality],
        eukaryotic: &[bool],
    ) -> Option<Self> {
        let hosts = offered
            .found
            .iter()
            .zip(scores.iter().zip(eukaryotic))
            .filter(|(_, (quality, eukaryotic))| Role::of(**quality, **eukaryotic) == Role::Host)
            .map(|(contigs, _)| judge.mapped(contigs))
            .collect::<Vec<_>>();
        let agreement = Agreement::of(
            &hosts,
            judge.depth.depths(),
            judge.depth.composition(),
            judge.lengths,
        )?;
        Some(Self {
            agreement,
            sums: HashMap::new(),
            joined: Vec::new(),
            taken: Vec::new(),
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
        rivals.sort_unstable_by_key(|(bases, rival, _)| (Reverse(*bases), *rival));
        let mine = self.agreement.sums(judge.mapped(piece));
        rivals.into_iter().find_map(|(_, rival, shared)| {
            let whole = self
                .sums
                .entry(rival)
                .or_insert_with(|| self.agreement.sums(judge.mapped(&offered.found[rival])));
            let rest = whole.without(&self.agreement.sums(judge.mapped(&shared)));
            self.agreement.agree(&mine, &rest).then_some(rival)
        })
    }

    pub(crate) fn record(&mut self, cluster: usize, at: usize, target: Option<usize>) {
        self.taken.push((cluster, at));
        if let Some(target) = target {
            self.joined.push((cluster, target));
        }
    }

    pub(crate) fn merge(
        &self,
        cluster_map: BTreeMap<usize, Vec<usize>>,
    ) -> BTreeMap<usize, Vec<usize>> {
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
        let mut groups = BTreeMap::<usize, Vec<usize>>::new();
        for (cluster, members) in cluster_map {
            let head = root(&mut parent, cluster);
            groups.entry(head).or_default().extend(members);
        }
        debug!(
            "Under-floor joins: {} pieces onto {} candidates, {} bins left",
            self.joined.len(),
            first.len(),
            groups.len()
        );
        groups
            .into_values()
            .enumerate()
            .map(|(id, mut members)| {
                members.sort_unstable();
                (id, members)
            })
            .collect()
    }
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
