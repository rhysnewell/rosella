use std::cmp::Reverse;
use std::collections::{BTreeMap, HashMap};

use log::debug;

use crate::quality::Quality;
use crate::recover::agree::{Agreement, Sums};
use crate::recover::fold::Role;
use crate::recover::join_report::Report;
use crate::recover::ladder::{Judge, Offered};

// Under the floor the smaller candidate is taken first, so a genome cut on depth noise leaves as
// pieces. A piece that agrees with the rest of a larger candidate holding it rejoins that one.
pub(crate) struct Joins<'a> {
    agreement: Agreement<'a>,
    sums: HashMap<usize, Sums>,
    whole: HashMap<usize, bool>,
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
            whole: HashMap::new(),
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
        rivals.sort_unstable_by_key(|(bases, rival, _)| (Reverse(*bases), *rival));
        let mapped = judge.mapped(piece);
        let mine = self.agreement.sums(mapped.iter().copied());
        if let Some(report) = &mut self.report {
            report.piece(&mapped);
        }
        let mut chosen = None;
        for (bases, rival, shared) in rivals {
            let found = judge.mapped(&offered.found[rival]);
            let whole = self
                .sums
                .entry(rival)
                .or_insert_with(|| self.agreement.sums(found.iter().copied()));
            let rest = whole.without(&self.agreement.sums(judge.mapped(&shared)));
            let agree = self.agreement.agree(&mine, &rest);
            let union = match shared.len() == piece.len() {
                true => found.clone(),
                false => {
                    let held = shared.iter().collect::<std::collections::HashSet<_>>();
                    let outside = piece.iter().filter(|at| !held.contains(at)).copied();
                    found
                        .iter()
                        .copied()
                        .chain(judge.mapped(&outside.collect::<Vec<_>>()))
                        .collect()
                }
            };
            let one = match (agree, shared.len() == piece.len()) {
                (false, _) => false,
                (true, true) => *self
                    .whole
                    .entry(rival)
                    .or_insert_with(|| self.agreement.holds_one(&union)),
                (true, false) => self.agreement.holds_one(&union),
            };
            let pick = chosen.is_none() && agree && one;
            if pick {
                chosen = Some(rival);
            }
            if let Some(report) = &mut self.report {
                let statistics = self.agreement.statistics(&mine, &rest);
                let (union_composition, union_depth) = self.agreement.scatter(&union);
                let (piece_composition, piece_depth) = self.agreement.scatter(&mapped);
                let rest_contigs = union
                    .iter()
                    .copied()
                    .filter(|contig| !mapped.contains(contig))
                    .collect::<Vec<_>>();
                let quality = |contigs: &[usize]| judge.quality.score(contigs);
                let (p, r, u) = (quality(&mapped), quality(&rest_contigs), quality(&union));
                let fields = [
                    size as f64,
                    bases as f64,
                    judge.bases(shared.iter().copied()) as f64,
                    statistics.depth.0,
                    statistics.depth.1 as f64,
                    statistics.composition.0,
                    statistics.composition.1 as f64,
                    union_composition,
                    union_depth,
                    piece_composition,
                    piece_depth,
                    self.agreement.bars.0,
                    self.agreement.bars.1,
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
