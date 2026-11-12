use std::collections::{BTreeMap, BTreeSet, HashMap};

use crate::clustering::clusterer::Partitioning;
use crate::recover::ladder::candidates;
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::disjoint_set::DisjointSet;

impl RecoverEngine {
    // Every pass offers its bins, but a pass takes a place from the kept one only beyond the
    // markers' noise, as the walk does. Picked bin by bin, noise alone trades whole genomes away.
    pub(super) fn contest_passes(
        &self,
        passes: &[Partitioning],
        kept: usize,
        contigs: &[usize],
    ) -> Partitioning {
        let found = candidates(passes);
        let mut cluster_map = BTreeMap::new();
        for place in places(&found) {
            let held = |pass: usize| {
                place
                    .iter()
                    .filter(|at| found[**at].1.contains(&pass))
                    .map(|at| &found[*at].0)
                    .collect::<Vec<_>>()
            };
            let mut best = held(kept);
            if place.iter().any(|at| found[*at].1.len() < passes.len()) {
                let mut worth = self.bins_worth(&best, contigs);
                for pass in (0..passes.len()).filter(|pass| *pass != kept) {
                    let next = held(pass);
                    let gain = self.bins_worth(&next, contigs) - worth;
                    if gain > self.bins_spread(&best, &next, contigs) {
                        worth += gain;
                        best = next;
                    }
                }
            }
            for members in best {
                cluster_map.insert(cluster_map.len(), members.clone());
            }
        }
        let binned = cluster_map.values().flatten().collect::<BTreeSet<_>>();
        let outliers = passes
            .iter()
            .flat_map(|pass| {
                pass.outliers
                    .iter()
                    .chain(pass.cluster_map.values().flatten())
            })
            .filter(|at| !binned.contains(at))
            .copied()
            .collect::<BTreeSet<_>>();
        let contested = Partitioning {
            cluster_map,
            outliers: outliers.into_iter().collect(),
            score: None,
            arm: passes[kept].arm,
            seed: passes[kept].seed,
        };
        if let Some(report) = &self.stages {
            let names = &self.tables.coverage.contig_names;
            for (at, pass) in passes.iter().enumerate() {
                report.partitioning(&format!("passes_offered{at}"), pass, contigs, names);
            }
            report.partitioning("passes_combined", &contested, contigs, names);
        }
        contested
    }
}

// Bins that share a contig sit in one place, so a place is where the passes disagree.
fn places(found: &[(Vec<usize>, Vec<usize>)]) -> Vec<Vec<usize>> {
    let mut joined = DisjointSet::new(found.len());
    let mut first = HashMap::new();
    for (at, (members, _)) in found.iter().enumerate() {
        for contig in members {
            let (one, other) = (
                joined.find(at),
                joined.find(*first.entry(*contig).or_insert(at)),
            );
            if one != other {
                joined.attach(one, other);
            }
        }
    }
    let mut places = BTreeMap::<usize, Vec<usize>>::new();
    for at in 0..found.len() {
        places.entry(joined.find(at)).or_default().push(at);
    }
    places.into_values().collect()
}
